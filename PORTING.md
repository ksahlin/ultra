# uLTRA — Rust rewrite

## Goal

Port uLTRA from Python to Rust. The priorities, in the author's order:

1. **Easier to install**
2. **More stable**
3. **Faster**
4. **Lower peak memory**
5. **Bugfixes**
6. **Slim the repository** (added during reconnaissance; the repo is 2.1 GB and 99.2% of that is one directory)

The Python in `uLTRA` and `modules/` is the **normative reference**: when Rust and Python disagree,
Python is right until a human decides otherwise. Goals 3, 4 and 5 will sometimes conflict with
byte-identity. When they do the divergence is a **numbered Finding** stating what differs, why, and
what it costs — never a silent change.

This document is the deliverable of the reconnaissance session. Every number in it was **read off an
actual run on this machine** (Darwin 25.5.0, arm64), never derived by reimplementing the tool's
logic. Where a number is a cost measured on a deliberately mismatched input, that is said in place.

### The reference environment these numbers come from

```
conda create -n ultra_ref -c conda-forge -c bioconda \
  python=3.12 pip parasail-python python-edlib pysam dill intervaltree gffutils minimap2
```

Resolved to python 3.12.14, parasail 1.3.4, pysam 0.24.1, dill 0.4.1, gffutils 0.14,
minimap2 2.31-r1302. `namfinder` is **not on bioconda for osx-arm64** and was built from source
(`v0.1.3`, cmake 4.2.3, no workaround needed).

---

## Part 1 — Installation, which is goal 1 and is measured

### The three published paths

Measured with `conda create --dry-run` per subdir, `CONDA_OVERRIDE_GLIBC=2.17` so the
cross-platform solves are meaningful, miniforge3 on Darwin 25.5.0 / arm64.

| Published instruction | linux-64 | linux-aarch64 | osx-64 | osx-arm64 |
| --- | --- | --- | --- | --- |
| README, bioconda: `ultra_bioinformatics` | resolves (0.1) | resolves (0.1) | resolves (0.1) | **fails** |
| `INSTALL.sh <dir>` | **broken** | **broken** | **broken** | **broken** |
| README, from source: `pip install` the six deps | wheels | *see below* | wheels | **parasail source build** |

**osx-arm64 is the only subdir where the bioconda recipe cannot be solved**, and the reason is a
single missing package:

```
nothing provides edlib needed by ultra_bioinformatics-0.0.3.3-pyh5e36f6f_0
└─ ultra_bioinformatics 0.1 would require
   └─ namfinder =* *, which does not exist (perhaps a missing channel).
```

`namfinder` 0.1.3 is built for `linux-64`, `linux-aarch64` and `osx-64`, and **not for `osx-arm64`**.
That one gap is what makes the recipe unsolvable on Apple Silicon. Everything else it needs —
minimap2 2.31, parasail-python 1.3.4, pysam 0.24.1 — is present on all four subdirs.

### `INSTALL.sh` is broken on every platform, and it is two lines

Measured by running the script's own namfinder steps verbatim:

```sh
cmake -B build -DCMAKE_C_FLAGS="-march=native" -DCMAKE_CXX_FLAGS="-march=native"
make -j -C build
mv namfinder $path          # <- fails
```

`make -C build` puts the binary at **`build/namfinder`**. There is no `namfinder` at the repository
root — confirmed by `find . -name namfinder -type f`, which returns exactly `./build/namfinder`. So
`mv namfinder $path` fails, and because the script opens with `set -e` the whole installation aborts
there. The same bug sits in the minimap2 block's sibling logic. `-march=native` additionally makes
the resulting binary non-portable, which matters for anyone building once and shipping to a cluster.

### PyPI wheel coverage, which is what `pip install` actually hits

| package | latest | wheels | aarch64 / arm64 wheels |
| --- | --- | --- | --- |
| `parasail` | 1.3.4 | 5 | **NONE** |
| `edlib` | 1.3.9.post1 | 38 | macOS `universal2` only — **no linux aarch64** |
| `pysam` | 0.24.1 | 42 | macOS arm64, manylinux aarch64, musllinux aarch64 — **full coverage** |
| `dill`, `intervaltree`, `gffutils` | — | 1 each | pure Python, `py3-none-any` |

> **Hypothesis tested: "pysam is the worst install dependency." It is not — it is the best-covered
> compiled dependency in the set.** `parasail` is the worst: no aarch64 wheel of any kind, so both
> Apple Silicon and ARM Linux fall back to a source build. This is the same conclusion the
> NGSpeciesID port reached about the same package.

### What this means for the port's distribution target

A Rust binary removes `parasail`, `edlib`, `pysam`, `dill`, `intervaltree` and `gffutils` from the
install surface completely. What it **cannot** remove is decided by the subprocess audit below:
`minimap2` and `namfinder` are the only two external programs the live code paths invoke.

---

## Part 2 — Where the time and memory actually go

### Corpora used

| corpus | reference | reads | note |
| --- | --- | --- | --- |
| smoke | `test/SIRV_genes.fasta` (7 seqs, 225 859 B) + `test/SIRV_genes_C_170612a.gtf` (339 lines) | `test/reads.fa` — **4 reads** | the repo's committed fixture |
| sirv-10k | same | 10 000 real ONT reads, 14 155 093 B | |
| sirv-100k | same | 100 000 real ONT reads, 141 954 278 B | the workhorse |
| gencode | — | — | GENCODE v47 human GTF, **1 805 576 649 B**, 4 105 490 lines |
| droso | `genomes/fruitfly.fa` (contigs named `2L`, `2R`, …) | `lrRNA-seq/droso/full_length_output_first_*.fq` | **needs a GTF that does not exist locally** — see below |

The only Drosophila annotation on hand is `annotations/dmel-all-r6.68.gff`, and it is unusable by
uLTRA as it stands: it is `##gff-version 3` with GFF3 attributes (`ID=`, `Parent=`) rather than
`gene_id`/`transcript_id`; its transcript feature is named **`mRNA`**, while uLTRA queries
`features_of_type('transcript')` and would therefore find **zero**; and it is 6.77 GB because the
genome is appended after `##FASTA` at line 31 790 958. Either convert with AGAT as the README
advises, or fetch the Ensembl BDGP6 GTF, whose contig naming already matches `fruitfly.fa`.

> **The committed fixture is four reads.** Rule 5 of this project exists because of corpora like
> this one. Every parameter in the tool is invisible to it.

### `index`

| corpus | wall | peak tree RSS | output |
| --- | --- | --- | --- |
| SIRV (339-line GTF) | **0.32 s** | **62.9 MB** | 608 KB total: `database.db` 245 760 B + 20 pickles, 368 KB |

### `align`, SIRV reference, real ONT reads

| reads | `--t` | wall | peak tree RSS | `reads.sam` |
| --- | --- | --- | --- | --- |
| 10 000 | 3 (default) | 10.95 s | 1 147.0 MB | 16 428 520 B |
| 100 000 | 1 | 253.76 s | 1 118.5 MB | |
| 100 000 | **3 (default)** | **92.49 s** | **1 914.8 MB** | 164 719 944 B |
| 100 000 | 8 | 44.58 s | **3 474.6 MB** | |
| 100 000 | 3, `--disable_mm2` | 58.63 s | 1 882.5 MB | |

### Where the 92.5 s goes — the tool's own stage timers, 100 000 reads, `--t 3`

| stage | time | share |
| --- | --- | --- |
| minimap2 (external subprocess) | 30.53 s | 33.0 % |
| processing reads | 2.15 s | 2.3 % |
| namfinder (external subprocess) | 1.32 s | 1.4 % |
| **uLTRA's own alignment loop (Python)** | **54.98 s** | **59.4 %** |
| selecting final best alignments | 3.29 s | 3.6 % |

**This bounds goal 3 and it should be stated before any Rust is written.** 34 % of default wall clock
is minimap2, which a port does not make faster. If uLTRA's own 55 s went to *zero*, the 100k run
would go from 92.5 s to about 37.5 s — a **2.5× ceiling** at the default settings. The honest targets
are therefore: a large speedup on the 59 % that is uLTRA's own code, and the `--disable_mm2` path
(58.6 s, of which 55 s is ours) as the number that shows the port's real work.

### Memory: it grows with `--t`, and the index is not why

1 118 MB at `--t 1` → 1 915 MB at `--t 3` → 3 475 MB at `--t 8`: about **337 MB per additional
worker**, on a reference whose entire index is **608 KB**. So on this corpus the per-worker cost is
Python process overhead and buffered data, not index data.

On a mammalian index it is both, because of the mechanism below, and that is
[issue #15](https://github.com/ksahlin/ultra/issues/15) — a user with 48 cores and 196 GB
who has to throttle `--t` to avoid OOM.

---

## Part 3 — The three hypotheses, tested

### Hypothesis 1: "gffutils' sqlite index dominates `index` memory and startup." — **Half right, and the wrong half is the memory half.**

Measured directly on GENCODE v47 human (1.8 GB GTF, 4 105 490 lines: 2 155 005 exon, 385 659
transcript, 78 724 gene), using uLTRA's own `create_db` argument set:

| operation | time | peak RSS | on-disk |
| --- | --- | --- | --- |
| `gffutils.create_db(..., disable_infer_genes=True, disable_infer_transcripts=True)` | **140.6 s** | **174 MB** | **3 719 MB** |
| `db.features_of_type('exon', order_by='seqid')` → 2 155 005 rows | **18.2 s** | 226 MB | — |

So sqlite dominates `index` **time** (140 s before uLTRA's own work begins) and **disk** (a 3.7 GB
database), and is **cheap in memory** (174 MB). Note this is the *fast* path: `--disable_infer` is
set. The default path infers genes and transcripts and is slower still — measured on Drosophila at
**240.3 s against 21.1 s**, for a semantically identical index. See *Finding 23*.

And uLTRA's actual use of gffutils is **three query patterns**:

```python
db.features_of_type('exon', order_by='seqid')
db.features_of_type('transcript', order_by='seqid')
db.children(transcript, featuretype='exon', order_by='start')
```

A direct GTF parser in Rust satisfies all three with a sort, and needs no database, no sqlite and no
3.7 GB file. This is the single largest, lowest-risk win in the whole port.

### Hypothesis 2: "pysam is the worst install dependency." — **Wrong.**

pysam has the most complete aarch64 wheel coverage of any compiled dependency here (see the table in
Part 1). Its *usage* is also shallow — reading minimap2's SAM, writing two SAM files from a template,
and `qualities_to_qualitystring`. The worst install dependency is **parasail**.

### Hypothesis 3: "dill implies the parallelism is pickling large objects between processes." — **Two separate things, and the stated form is wrong.**

`dill` is used for the **on-disk index**, not for inter-process communication.
`help_functions.pickle_dump/pickle_load` write and read the 20 index pickles with dill.

Why dill and not stdlib pickle: `create_augmented_gene.py` builds
`defaultdict(lambda: array("L"))`, and a lambda is not picklable by the standard library. Proven
rather than assumed — loading `parts_to_segments.pickle` with **stdlib `pickle`** succeeds but
**imports `dill` as a side effect of the load**, because the stream itself references dill's
reconstruction helpers:

```
loaded with stdlib pickle, dill imported as side effect: True
default_factory: <function get_canonical_segments.<locals>.<lambda> at 0x101216980>
```

So dill is a genuine runtime dependency, and it is caused by four lambdas. Replacing them with
module-level named functions would drop dill from the Python tool entirely — a worthwhile upstream
commit independent of the port.

**But there *is* heavy IPC pickling**, just not dill's: `pc.py` passes batches of 1 000 reads with
all of their seeds through `mp.Manager().Queue()`, which pickles with the standard library. And each
worker calls `align.import_data(args)` at start-up, so **the entire index is loaded once per `--t`
process**:

```python
def align_single(process_id, input_queue, output_sam_buffer, ...):
    auxillary_data = import_data(args)      # all 13 pickles, per worker
```

That is the mechanism behind issue #15, and it is the thing a Rust port removes for free: one
immutable index, shared by reference across threads.

---

## Part 4 — Dependency and subprocess surface

### Python packages

| dependency | used for | verdict |
| --- | --- | --- |
| `parasail` | 4 call sites in `help_functions.py`: `sg_trace_scan_16`/`_32`, `ssw`, `sw_trace_scan_16`, with `matrix_create("ACGT", …)` | **Rust crate or vendored C.** `libparasail-sys` is exact by construction; the NGSpeciesID port also proved a hand-written `parasail.rs` can be exact. Decide by oracle, not by preference |
| `edlib` | 2 call sites: `classify_read_with_mams.py` and `help_functions.py`, `task="path"` | **Vendor edlib's single `.cpp` with `cc`**, or reimplement against a recorded oracle. The NGSpeciesID port found `edlib_rs` and `rsedlib` both unusable |
| `pysam` | read minimap2's SAM; write `indexed.sam`/`unindexed.sam` from a template; `qualities_to_qualitystring` | **Replace.** `noodles-sam`, or a hand-rolled writer — the header uLTRA emits is 7 `@SQ` lines |
| `dill` | the 20 on-disk index pickles | **Gone.** The Rust index will not be a Python pickle — see *Finding 3* |
| `gffutils` | three query patterns, over a 3.7 GB sqlite | **Gone.** Direct GTF parser |
| `intervaltree` | one use, `prefilter_genomic_reads.py`, one tree per chromosome of uLTRA-indexed regions | **Replace.** `coitrees` / `rust-lapper`, or a sorted vec + binary search |

### External programs — audited, and two of the four documented ones are dead code

| program | invoked from | live? |
| --- | --- | --- |
| `namfinder` | called from **`uLTRA:364`**, unconditionally, via `seed_wrapper.find_nams_namfinder`, which uses **`os.system`** with a shell pipeline: `namfinder … -S ref reads 2> err \| gzip -1 --stdout > seeds.txt.gz` | **live, and it is the *only* seed finder** — see *The seeding decision* |
| `minimap2` | `prefilter_genomic_reads.align_with_minimap2`, `subprocess.check_call` | **live — must stay a subprocess** |
| `mummer` | `seed_wrapper.find_mems_mummer` | **dead** — called only from `find_mems_slamem` |
| `slaMEM` | `seed_wrapper.find_mems_slamem` | **dead** — nothing calls it |

`get_mem_records` — the only consumer of the mummer/slaMEM output — is also dead: nothing calls it,
and it would raise `NameError` if anything did, because the `mem` namedtuple it uses is commented out
at the top of `seed_wrapper.py`. The two live call sites use `read_seeds` instead.

**So the port's external-tool constraint is exactly two binaries, both of which must be on `PATH`.**
That is what bounds the "portable prebuilt binary" story: the uLTRA binary itself can be static, but
as things stand it will still shell out to `minimap2` and `namfinder`.

### The seeding decision — open

There is **no flag selecting a seed finder**. `uLTRA:364` calls namfinder unconditionally, and the
MEM-based alternatives are the dead code of Finding 7. So strobemer/NAM seeding is not one option
among several; it is what the tool does.

That matters more than it looks, because **namfinder is the sole reason the bioconda recipe fails on
osx-arm64** (Finding 2). Three ways out, and they are not equivalent:

| option | install effect | byte-identity | scope |
| --- | --- | --- | --- |
| **A. keep namfinder as a subprocess** | osx-arm64 stays broken until the feedstock gains a build | achievable — same binary, same seeds | a port |
| **B. compute strobemers inside the Rust binary** | **removes an external tool and fixes osx-arm64 outright**; only minimap2 left | achievable *if* the NAM output is reproduced exactly — an oracle problem of the same shape as spoa in the NGSpeciesID port | a port, plus one hard oracle |
| **C. revive MEM seeding** | still needs mummer or slaMEM, or a new MEM finder | **impossible** — different seeds give different alignments | a research change, not a port |

**Decision taken: B, and implemented by LINKING rather than reimplementing.** Strobemers are computed
inside the Rust binary and namfinder is no longer a runtime dependency — but the port does not
reimplement randstrobe seeding. It vendors namfinder's own sources and calls `run_strobealign`, the
very function namfinder's `main` calls. See *Finding 26*: this is exact by construction, so the NAM
oracle became a regression check rather than a specification, and fallback A was never needed.

C is *Deferred improvements*. MEMs may well be more sensitive than NAMs, but changing the seeds
changes every output file, so C cannot be evaluated until the port is exact enough to serve as the
control. Port first, then measure MEMs against a known-good baseline. Reviving it would also mean
reviving `get_mem_records`, which is dead *and* broken (Finding 7).

This is the treatment spoa got in the NGSpeciesID port, and for the same reason: *"a vendored library
built at compile time is strictly better than either alternative. It needs a C++ toolchain and CMake
to build, and nothing at all to run."* Here it does not even need CMake.

Filed upstream: **[ksahlin/namfinder#1](https://github.com/ksahlin/namfinder/issues/1)** asks for
`osx-arm64` in the bioconda recipe. That fixes the *reference's* install on Apple Silicon
independently of the port, and it stays worth doing even under B, because the reference environment
the harness runs against still needs a namfinder binary.

---

## Part 5 — The equivalence contract

### SAM is eligible as a byte-identity target. Measured.

The author flagged `@PG` as a risk, because NGSpeciesID's medaka BAM hit exactly that (its Finding
24). **uLTRA's own output does not have the problem:**

| file | `@HD` | `@PG` | verdict |
| --- | --- | --- | --- |
| **`reads.sam`** (the deliverable) | absent | **absent — 0 lines** | **byte-identity target** |
| `minimap2.sam`, `indexed.sam`, `unindexed.sam` | present | present, carries `VN:` and the full `CL:` command line | temp files; see below |

`reads.sam` is written by
`pysam.AlignmentFile(name, "w", reference_names=..., reference_lengths=...)`, which emits `@SQ` lines
and nothing else. No program record, no version, no paths.

The three minimap2-derived files **do** carry `@PG` with minimap2's version and the verbatim command
line — including whatever paths the user passed — but all three are **deleted unless
`--keep_temporary_files`**. The decision is therefore cheap: hash `reads.sam` unconditionally, and
under `--keep_temporary_files` hash the minimap2 files with the `@PG` line scrubbed, exactly as
NGSpeciesID's `NOT_CONTRACT_RE` handles racon's stderr.

### Determinism: `reads.sam` is reproducible; the index is not

Two full `pipeline` runs with `--keep_temporary_files`, all 34 output files hashed:

| condition | files differing |
| --- | --- |
| default environment | **17 of 34** |
| `PYTHONHASHSEED=0` | **3 of 34** |

`reads.sam` is **identical in both conditions** — it is deterministic even under hash randomisation.

The 14 files that `PYTHONHASHSEED=0` fixes are the index pickles and `refs_sequences.fa`: they
serialise sets and dicts whose iteration order is hash-seed dependent. The 3 that remain are fully
explained:

| file | cause | content |
| --- | --- | --- |
| `reads_tmp.fa.gz`, `seeds.txt.gz` | gzip **mtime**, header bytes 5–8 (`44bf ae6a` vs `45bf ae6a`, one second apart) | decompressed streams **byte-identical** |
| `minimap2_errors.1` | minimap2's own timing lines | — |

**So the harness must set `PYTHONHASHSEED=0` for every reference run**, exactly as NGSpeciesID's
does, and must compare the two `.gz` files decompressed.

### The index cannot be a byte-identity target, and that is structural

The Python index is 20 **Python pickles** plus a **gffutils sqlite database**. A Rust port will write
neither. This is the one place where uLTRA differs fundamentally from the NGSpeciesID port, whose
every output was text.

The contract is therefore split:

- **`reads.sam` — byte-identical.** This is the acceptance criterion.
- **the index — semantically identical**, verified by a stage oracle: `bench/dump_reference.py`
  loads the Python pickles and emits a canonical, sorted, text rendering of each of the 14
  structures; the Rust binary emits the same rendering; the two are diffed. This is the
  "stage oracles are required, not optional" rule, and here it is not optional at all.

---

## Part 6 — Findings in the reference

Numbered as they were found. Each is a decision to make, not a bug to fix silently.

### Finding 1 — `INSTALL.sh` cannot succeed on any platform
`mv namfinder $path` runs at the repository root; cmake writes `build/namfinder`; `set -e` aborts.
Measured by building namfinder v0.1.3 from source. Also builds with `-march=native`, which makes the
binary non-portable. **Cost to fix: two lines.** Belongs on a `fix/installation` branch that does not
wait for the port.

### Finding 2 — the bioconda recipe is unsolvable on osx-arm64, because namfinder is not built for it
`namfinder` 0.1.3 exists for linux-64, linux-aarch64 and osx-64. Adding an osx-arm64 build to the
namfinder feedstock fixes uLTRA's install on Apple Silicon without touching uLTRA. **This is the
highest-leverage single action for goal 1 and it is in a different repository.**

### Finding 3 — the index format cannot be preserved, so `index` output is not a byte-identity target
20 Python pickles + a gffutils sqlite database. Documented above; drives the stage-oracle design.

### Finding 4 — the index is loaded once per `--t` worker
`align_single` calls `import_data(args)` per process. Memory scales as `--t × index`. This is
[issue #15](https://github.com/ksahlin/ultra/issues/15). The Rust port shares one immutable index
across threads and the problem disappears — this is the strongest argument for the port after
installation.

### Finding 5 — `buffer_write_cnt` is never incremented, so the mid-run SAM flush never fires
`pc.py` sets `buffer_write_cnt = 0` at line 34, tests `if buffer_write_cnt >= 50000:` at line 51 and
resets it at line 53. **Nothing ever increments it.** The periodic drain of `output_sam_buffer` is
therefore dead, and every SAM record produced by every worker accumulates in the manager queue until
`file_IO` finishes reading input. Contributes to the peak RSS measured above. Faithful-port question:
reproduce the dead code, or fix it and take the divergence as a documented Finding.

### Finding 6 — `zip()` over reads and seeds silently truncates
`pc.file_IO` does
`zip(help_functions.readfq(open(reads)), seed_wrapper.read_seeds(seeds))`, and `read_seeds` yields
only when both a forward and a reverse record have been seen (`if curr_acc and curr_acc_rev`). `zip`
stops at the shorter iterator, so reads past the first mismatch are dropped **with no error**. This
is the mechanism behind [issue #3](https://github.com/ksahlin/ultra/issues/3), which the author
already diagnosed as hard to fix in the current parser. In Rust it is an error return, not a silent
truncation — and it is a behaviour change that needs its own Finding.

### Finding 7 — `mummer`, `slaMEM` and `get_mem_records` are dead code
Documented in Part 4. `get_mem_records` would `NameError` on `mem(...)` if reached. Deleting them
removes two tools from the README's dependency list at zero cost.

### Finding 8 — namfinder is invoked through `os.system` with an unquoted shell string
`seed_wrapper.py:74`. Paths are interpolated into a shell command with no quoting, so a space or a
shell metacharacter in `--index`, the reference path or the read path breaks or misbehaves. The Rust
port uses an argument vector and a pipe, which removes the class.

### Finding 9 — `.travis.yml` is dead CI pinned to Python 3.6–3.8 and tools the code no longer uses
Travis CI is gone; the file installs `slaMEM`, `StrobeMap` and `mummer`, none of which the live code
paths call. Replace, do not migrate.

### Finding 10 — the committed test corpus is four reads
`test/reads.fa`. Blind to essentially every parameter. See *Corpora*.

### Finding 11 — `--ont` and `--isoseq` advertise a flag that does not exist
Both help strings say they set `--min_seed` (18 and 20 respectively). **There is no `--min_seed`
flag**: the string appears only in those two help texts, and passing it is an argparse error —
`uLTRA: error: unrecognized arguments: --min_seed`. Both also claim to set
`--alignment_threshold 0.5`, which they do not (0.5 is merely the default). And neither documents
`mm2_ksize`, which `--ont` *does* change, from 15 to 14. The help text is wrong in three directions
at once; the port must decide whether to reproduce it verbatim or correct it as a Finding.

### Finding 12 — the presets silently override an explicit `--s`
Measured by reading the namfinder command line the tool prints:

| invocation | namfinder receives |
| --- | --- |
| `--s 12` | `-k 12 -s 12` |
| `--ont --s 12` | **`-k 9 -s 9`** |
| `--isoseq --s 12` | **`-k 10 -s 10`** |

`args.s` is overwritten after `parse_args()`, so a user who combines a preset with an explicit `--s`
gets the preset with no warning. Note also that `--isoseq` sets `s = 10`, which is already the
default, so that half of the preset is a no-op. Same class as NGSpeciesID's Finding 14, and each
preset needs its own equivalence case for exactly this reason.

### Finding 13 — `index` crashes on a standard Ensembl GTF, and the smoke corpus cannot see it

**The first corpus with a reference other than SIRV broke the reference implementation.** Indexing
*Drosophila melanogaster* (Ensembl BDGP6.46 release 113, 167 MB, 24 278 genes / 41 610 transcripts /
196 625 exons) against `genomes/fruitfly.fa`:

```
File "modules/create_augmented_gene.py", line 440, in create_graph_from_exon_parts
    assert active_start <= exon.start - 1
AssertionError
```

Exit 1 after 15.7 s. Root cause, measured rather than inferred — `create_augmented_gene.py:311`
iterates

```python
for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
```

and `order_by='seqid'` sorts by **contig only**. It says nothing about order *within* a contig, so
the rows come back in the database's natural order, and the loop body assumes ascending `start`.
Counting the exons gffutils actually returns:

| corpus | exons returned | returned out of ascending-start order |
| --- | --- | --- |
| SIRV (`test/SIRV_genes_C_170612a.gtf`) | 339 | **0** |
| Drosophila (Ensembl BDGP6.46) | 196 625 | **99 961** |

SIRV's annotation happens to be coordinate-sorted within each contig, so the assumption has never
been violated by anything in the repository. Ensembl's is grouped by gene, so it is violated 99 961
times.

**The fix is one word**: `order_by='seqid,start'`. Measured consequences:

| | |
| --- | --- |
| Drosophila index | **exit 1 → exit 0**, 9.25 s, 785.8 MB peak RSS, 640 MB index |
| SIRV index files changed | **2 of 21** — `parts_to_segments.pickle`, `gene_to_small_segments.pickle` |
| **SIRV `reads.sam`, 10 000 real ONT reads** | **byte-identical** |

So it converts a crash into a working index and does not change alignments. It is **not** a no-op on
the index itself, though, so it must not land between recording goldens and verifying against them —
and the byte-identity check above is so far one corpus, not five.

The sibling loop at line 481 uses the same `order_by='seqid'` but is order-independent (it
accumulates a max and a min, and its inner exon loop already passes `order_by='start'`), so the bug
is confined to line 311.

This is very likely related to [issue #20](https://github.com/ksahlin/ultra/issues/20), which is a
GENCODE-based index report, though that user's visible error comes from a downstream `gffcompare`
step and the connection is not yet proven.

Belongs on `fix/installation` — or rather a `fix/exon-ordering` branch, since it is a correctness fix
and not an installation one, and it changes index bytes.

### Finding 14 — `gene_to_small_segments` is built from a leaked loop variable, and is silently wrong on any real annotation

Found by the case matrix, not by reading: `idx-small-exon-50` recorded **exit 1** on SIRV.

```
File "modules/create_augmented_gene.py", line 251, in get_canonical_segments
    add_items(gene_to_small_segments[gene_id], chr_id, e_start, e_stop)
UnboundLocalError: cannot access local variable 'gene_id' where it is not associated with a value
```

`get_canonical_segments` tests `<= small_segment_threshold` in **three** places. Two of them are
written correctly:

```python
if p2 - p1 <= small_segment_threshold:
    for gene_id in active_gene_ids:                       # <- binds gene_id
        add_items(gene_to_small_segments[gene_id], chr_id, p1, p2)
```

The third omits the loop and uses a bare `gene_id`:

```python
if e_stop - e_start <= small_segment_threshold:
    add_items(gene_to_small_segments[gene_id], chr_id, e_start, e_stop)   # <- line 251
```

So it reads whatever `gene_id` Python left bound from an **earlier, unrelated** iteration of one of
the other two loops. `--small_exon_threshold 50` crashes only because at that threshold the earlier
loops do not run first, leaving the name unbound.

**The crash is the harmless symptom. The silent version is the bug.** Instrumenting the site to print
the leaked `gene_id` alongside the `active_gene_ids` actually in scope, at the **default**
`--small_exon_threshold 200`:

| corpus | site reached | leaked gene not in `active_gene_ids` | more than one active gene |
| --- | --- | --- | --- |
| SIRV | 22 | **0** | **0** |
| Drosophila (Ensembl BDGP6.46) | 1 048 | **6** | **835** |

On SIRV every part has exactly one active gene, so the leaked value is always the only correct value
and the bug is invisible. On a real annotation it is wrong twice over:

- **6 segments are filed under a gene that is not active there at all** — and one of the leaked
  values is not even a gene id but a transcript id, `FBtr0309646_df_nrg`:
  ```
  gene_id='FBgn0035170' in_active=False active=['FBgn0035171', 'FBti0020010']
  gene_id='FBtr0309646_df_nrg' in_active=False active=['FBgn0259163']
  ```
- **835 segments reach a part with several active genes** and are filed under one of them instead of
  all, so `gene_to_small_segments` is incomplete for the rest.

`gene_to_small_segments` is one of the 14 pickled index structures and feeds the small-exon handling
that is uLTRA's headline claim, so this is not cosmetic.

The fix is one line — wrap site three in `for gene_id in active_gene_ids:` like its two siblings. It
**will** change index bytes and may change alignments, so it is a behaviour change needing its own
branch, its own goldens and a measured before/after, not a drive-by.

For the port: the Rust code cannot reproduce "whatever was left in a local variable" and should not
try. This is a *deliberate divergence* — implement the evident intent (all active genes), and record
the diff against the reference on every corpus as the cost.

### Finding 15 — `--alignment_threshold` is `type=int` with a float default, so it cannot be set to the value its own help describes

Found by the case matrix: `aln-aln-thresh-03` recorded **exit 2**.

```
uLTRA align: error: argument --alignment_threshold: invalid int value: '0.3'
```

Both declarations say:

```python
add_argument('--alignment_threshold', type=int, default=0.5, help='... Default val (0.5) sets that a
             score higher than 2*0.5*read_length would be considered an alignment ...')
```

`argparse` does not run `type` over a default that is not a string, so the **float default of 0.5 is
used as-is and works**, while **anything the user passes must parse as an `int`**. The help text
describes a fraction and gives 0.5 as the example, so the one value the documentation points at is
the one value that cannot be typed on the command line. `--alignment_threshold 1` is accepted and
means something entirely different.

Exit code 2 is argparse's, so this is a clean CLI-contract case rather than a crash. `--max_loc` is
declared `type=float, default=5` — harmless in the other direction, and worth keeping an eye on for
the same reason.

The port must decide: reproduce the `int` restriction exactly (and keep a flag nobody can use
correctly), or make it `float` and take the divergence. Recommend `float`, as its own commit, with
the CLI golden for `aln-aln-thresh-03` changing from exit 2 to exit 0 as the documented cost.

### Finding 16 — an invalid `--thinning` prints an error and exits **0**

```
$ uLTRA align --thinning 3 ref.fa reads.fa out/
Invalid thinning level. Choose 0, 1 or 2.
$ echo $?
0
```

`uLTRA:551` does `print(...)` then a bare `sys.exit()`, and `sys.exit()` with no argument exits **0**.
So a rejected input is indistinguishable from a successful run to any caller that checks the exit
status — a Snakemake rule, a nextflow process, `set -e`. Same for `--thinning -1`.

Contrast with the argparse failures, which are correct: `--t abc`, `--ont --isoseq` and an unknown
flag all exit 2.

**And it happens in a second place**: a missing `--index` folder is detected, reported clearly, and
also exits 0 (see Finding 18). So two distinct "user got it wrong" paths both report success.

One line each: `sys.exit(1)`. It is a behaviour change, so both get their own CLI goldens.

### Finding 17 — failures after argument parsing leave an empty output folder behind

`parse_args()` is at `uLTRA:528` and `help_functions.mkdir_p(args.outfolder)` at `uLTRA:533`, so the
side effect splits cleanly along that line. Measured across all 25 CLI cases:

| failure class | exit | outfolder |
| --- | --- | --- |
| argparse (`--t abc`, `--ont --isoseq`, unknown flag, bad subcommand, missing positionals) | 2 | **not created** |
| the `--thinning` range check | 0 | **created, empty** |
| missing index folder | 0 | **created, empty** |
| runtime failures (missing reference / reads / gtf) | 1 | **created, empty** |

So a user who mistypes a flag gets nothing on disk, and a user who mistypes a *path* gets a stray
empty directory. It is observable, so the CLI goldens record `created=` and `entries=` next to the
exit code.

### Finding 18 — `align` against a nonexistent reference blames the wrong file

```
$ uLTRA align /nonexistent.fa reads.fa out/
FileNotFoundError: ... '/tmp/out/ref_part_sequences.pickle'
```

Exit 1 with a traceback. The reference path is never checked, so the first thing that actually fails
is the index load, and the error names an index pickle the user has never heard of rather than the
reference they mistyped. The same shape of error appears for a missing index folder, which makes the
two indistinguishable.

Cheap to improve in the port (check inputs exist, name the one that is missing), and it is exactly
the kind of thing the *more stable* goal is about — but it changes stderr, so it is a Finding with a
CLI golden, not a silent fix.

Note the contrast with a missing **index folder**, which the tool does handle properly — it is
detected, named, and explained:

```
The index folder specified for alignment is not found. You specified:  /nonexistent_index_dir
Build  the index to this folder, or specify another forder where the index has been built.
```

…and then exits **0**, which is Finding 16 again in a second place. So the machinery for a good error
message already exists; it is the reference path that is unchecked, and the exit code that is wrong.
(Two typos in that message, `forder` and the double space, are contract until deliberately changed.)

### Finding 19 — argparse's help text depends on `$COLUMNS`, so it is machine-dependent output

`argparse` wraps usage and help to the terminal width, which it takes from `$COLUMNS` (falling back
to 80 when it cannot tell). Measured on `uLTRA --help`:

| `COLUMNS` | bytes | longest line |
| --- | --- | --- |
| unset | 761 | 78 |
| 80 | 761 | 78 |
| 100 | 713 | 97 |
| 200 | 689 | 110 |

So the *same* reference, on the *same* input, emits different bytes depending on an environment
variable the user probably does not know is set. Any CLI golden recorded without pinning it is a
false failure waiting for the first developer with a wide terminal.

`equivalence.sh` now sets `COLUMNS=80` on both the reference and the port. Re-recording all 25 CLI
goldens with the pin in place was a byte-for-byte no-op, which confirms 80 is what they were
originally captured at.

For the port this is a small liberation: the Rust binary emits **fixed** strings extracted from the
goldens, so it is correct at `COLUMNS=80` and does not reflow. Matching argparse's reflow at other
widths is not attempted and is not contract — if it ever needs to be, this Finding is where the
decision gets recorded.

### Finding 20 — the index keys are platform-width dependent, and that fixes the port's encoding

Every segment, part, exon and flank is keyed by `array("L", [chr_id, start, stop]).tobytes()` and
read back with `struct.unpack("LLL", key)`. Both use the **native** `unsigned long`, so:

| platform | `array("L").itemsize` | key length |
| --- | --- | --- |
| Linux / macOS, 64-bit (LP64) | 8 | 24 bytes |
| Windows 64-bit (LLP64), and any 32-bit | 4 | 12 bytes |

The two halves agree with each other, so nothing breaks *within* a platform — but an index built on
one cannot be read on another. In practice uLTRA is Unix-only, so this has never bitten anyone.

For the port this is good news rather than bad: **Rust `u64` little-endian is byte-identical to
`array("L")` on every platform uLTRA actually runs on**, so the port writes `u64` LE explicitly and
is both compatible and no longer platform-dependent. That is a silent improvement, and it is
recorded here rather than left to be rediscovered.

### Finding 21 — two index structures vary run-to-run, and it provably does not matter

`parts_to_segments` and `gene_to_small_segments` are the only `array("L")` structures, and their
contents are in **set-iteration order**, so they change with `PYTHONHASHSEED`. Measured on SIRV,
comparing two indexes built at seeds 0 and 12345:

| | differ |
| --- | --- |
| raw pickles | **14 of 20** |
| oracle rendering, arrays in stored order | 2 of 20 — exactly these two |
| oracle rendering, arrays sorted | **0 of 20** |

> **THIS FINDING WAS WRONG AND IS CORRECTED BY FINDING 27.** It is left here, with the mistake
> intact, because how it was wrong is more useful than the conclusion was.

The claim was that the order is noise, on two grounds:

- **Both consumers de-duplicate.** `gene_to_small_segments` is read into a `set(...)` at
  `classify_read_with_mams.py:266`. `parts_to_segments` feeds `segment_hit_locations`, which is
  `list(set(segment_hit_locations))` and then sorted at `:244-245`.
- **`reads.sam` is byte-identical across hash seeds** while 14 of the 20 pickles are not.

Both statements are true. **Neither supports the conclusion.**

De-duplicating into a `set` does not remove order: it replaces the input order with the *set's*
iteration order, which still depends on the insertion sequence. And the sort at `:245` has key
`x[1]` only, so entries tied on `x[1]` keep whatever order the set produced. Python's sort is stable,
which is exactly what lets the upstream order survive.

The second was a measurement taken where it could not fail: on the 4-read smoke corpus, and on runs
that reused one pre-built index, so the array order never actually varied between the things being
compared. Measured properly — two indexes built at different hash seeds, 10 000 real reads —
`reads.sam` **differs**. See Finding 27.

`bench/dump_reference.py` still sorts arrays, and that is now a deliberate weakening of the oracle
rather than a free canonicalisation: it means the stage-index comparison cannot see an order
difference that does reach output. Finding 27's fix makes the reference's order deterministic, which
is what restores the guarantee.

### Finding 22 — DIVERGENCE TAKEN: `gene_to_small_segments` is built for every active gene

The port does **not** reproduce Finding 14's leaked loop variable. It implements the evident intent:

```rust
// reference: add_items(gene_to_small_segments[gene_id], ...)   <- gene_id is whatever
//            was left bound by an earlier, unrelated loop
for gene_id in &active_gene_ids {
    add_items(gene_to_small_segments.entry(gene_id), chr_id, e_start, e_stop);
}
```

matching the two sibling call sites that were written correctly. Decision taken by the author on
2026-09-20: *implement the evident intent, document the divergence.*

**What it costs.** The port's `gene_to_small_segments` will differ from the reference's wherever a
part has more than one active gene, or where the leaked value was simply wrong. Measured at the
default `--small_exon_threshold`:

| corpus | site reached | reference files under a gene that is not active | reference reaches a part with >1 active gene |
| --- | --- | --- | --- |
| SIRV | 22 | 0 | 0 |
| Drosophila | 1 048 | 6 | 835 |

So on SIRV the divergence is expected to be **empty** — every part there has exactly one active gene,
which makes SIRV useless for confirming the fix works and essential for confirming it breaks nothing.
On Drosophila the two will differ on up to 841 entries, and the port's is the correct one.

**How it is verified rather than asserted.** `equivalence.sh stage index` compares the oracle
rendering structure by structure, so this divergence is confined to one named file
(`gene_to_small_segments.txt`) and every other structure must still match exactly. A divergence that
leaks into a second structure is a bug, not a decision. The expected-diff list lives in
`bench/stage_diffs.tsv` so it cannot be widened silently.

**Measured, now that the port exists.** Across all six corpora, 120 structures:
**119 match, 1 diverged, 0 failed, 0 stale.** The single divergence is this one, on droso-20k, and
its shape is exactly what the diagnosis predicted:

| | |
| --- | --- |
| SIRV `gene_to_small_segments` | **matches** — as it must; every SIRV part has one active gene |
| droso genes in the reference / the port | 16 640 / **16 642** |
| shared genes where the port is a superset | 16 635, **strictly larger in 836** |
| shared genes where the port has fewer entries | **5 genes, 6 entries** |

The 836 are the multi-gene parts the reference filed under only one gene (it measured 835 sites; one
gene is reached twice). The 6 entries the port "loses" are exactly the 6 the reference filed under a
gene that was **not active there at all** — the port is right to omit them, and the count matching
the independent instrumentation is what confirms the two measurements are describing one bug.

### Finding 23 — the default `index` path spends 219 of its 240 seconds producing an identical index

`--disable_infer` tells gffutils not to infer `gene` and `transcript` features. The README calls it a
speed-up "if you have the gene feature and transcript feature in your GTF file". Measured on
Drosophila (Ensembl BDGP6.46, 167 MB, which *does* carry 24 278 `gene` and 41 610 `transcript` rows)
against `genomes/fruitfly.fa`:

| path | wall | peak tree RSS | `database.db` | pickles |
| --- | --- | --- | --- | --- |
| default | **240.3 s** | 863.1 MB | 399 MB | 271 MB |
| `--disable_infer` | **21.1 s** | 827.0 MB | 388 MB | 270 MB |

**11.4× faster — and the index is semantically identical**: the stage oracle renders
**0 of 20 structures differently** between the two. The 11 MB of extra `database.db` is the inferred
features themselves, which nothing downstream reads.

So on any annotation that already declares its genes and transcripts — GENCODE, Ensembl, essentially
every modern GTF — the default path burns 219 seconds building rows that make no difference. The
flag is documented as a speed-up and is really closer to "do not waste four minutes".

Note the converse still holds: for a GTF that genuinely lacks `gene`/`transcript` rows, inference is
required and `--disable_infer` would produce a *different* and wrong index. So the port cannot simply
hardwire it. What it can do is **detect** the features and skip inference when they are present,
which is the same decision the user is currently asked to make by hand — and get it right by
construction.

The port replaces gffutils entirely (Finding 1 of the hypotheses, Part 3), so neither path survives
as such. This Finding exists to record what the target is: a GTF parse, not 240 seconds and a 399 MB
database.

### Finding 24 — a dedup guard compares a tuple to a bytes-keyed dict, and is dead

In `get_canonical_segments`, the "extend forwards" branch guards with

```python
if (chr_id, p1, pos_tuples[i+k][1]) not in segment_id_to_choordinates:
```

while its "extend backwards" sibling twenty lines earlier correctly uses the bytes key:

```python
if segment_name not in segment_id_to_choordinates:
```

`segment_id_to_choordinates` is keyed by `array("L",...).tobytes()`, so a **tuple** can never be in
it and the guard is always true. It is dead code.

I predicted from reading that this would produce duplicate entries, and **measured that it does
not**: zero duplicate triples in `parts_to_segments` and `gene_to_small_segments` on both SIRV
(275 / 221 triples) and Drosophila (102 065 / 66 790). Each `i` adds at most one forward segment
before `break`, and distinct `i` give distinct `p1`, so nothing collides in practice.

So it is a latent defect with no observable effect. The port reproduces the *behaviour* — always add
— rather than the *intent*, because implementing a working dedup here could drop a segment the
reference keeps, which would be a real divergence in exchange for nothing.

Worth recording as a method note too: this is the second time in this port that reading the code
predicted something the measurement contradicted.

### Finding 25 — `--disable_infer` on a GTF without transcript lines silently empties the annotation

`--disable_infer` tells gffutils not to infer `gene` and `transcript` features. **The repository's own
test GTF has neither** — `test/SIRV_genes_C_170612a.gtf` is 339 `exon` lines and nothing else.

Running the reference with `--disable_infer` on it:

```
$ uLTRA index --disable_infer test/SIRV_genes.fasta test/SIRV_genes_C_170612a.gtf out/
$ echo $?
0
```

and the resulting index has **every splice structure empty**:

| structure | default | `--disable_infer` |
| --- | --- | --- |
| `splices_to_transcripts` | 7 | **0** |
| `transcripts_to_splices` | 7 | **0** |
| `all_splice_pairs_annotations` | 7 | **0** |
| `all_splice_sites_annotations` | 7 | **0** |
| `max_intron_chr` | 7 | **0** |

Exit 0, no warning. uLTRA is an *annotation-guided* aligner, so this silently removes the guidance
that is the entire point of the tool — and Finding 23 actively encourages users towards this flag by
making it 11.4× faster.

The two findings together are the trap: the flag is a large speed-up on annotations that declare
their transcripts and a silent correctness disaster on annotations that do not, and nothing tells the
user which they have.

**The port does not offer the choice.** It always infers what is missing: a `transcript_id` seen on
an exon but never declared gets a transcript, one that was declared is left alone. There is no flag,
because there is no decision the user is better placed to make. `--disable_infer` is still accepted
and ignored.

### Finding 26 — namfinder already builds as a library, so the seeding stage is exact by construction

The seeding decision assumed the choice was *reimplement randstrobes* or *keep the subprocess*. It is
not. namfinder's own CMake already separates a static library from its CLI:

```cmake
add_library(salib STATIC ... src/randstrobes.cpp src/nam.cpp src/output.cpp ...)
add_executable(namfinder src/main.cpp)
target_link_libraries(namfinder PUBLIC salib)
```

and `main()` is a four-line wrapper around `run_strobealign(argc, argv)`. Everything that decides the
output — argument parsing, index construction, NAM finding, and `output_nams`, which emits the
`> read` / `> read Reverse` format uLTRA parses — is in the library.

So the port vendors `src/` and `ext/` (**1.1 MB, 51 files, MIT, same author**) and calls
`run_strobealign` through a six-line C ABI shim. `rust/build.rs` compiles exactly upstream's `salib`
source list with `cc`, plus `main.cpp` with `-Dmain=namfinder_cli_main` so its entry point does not
collide with Rust's. **No CMake is needed** — upstream's two generated headers carry one `#define`
each and `build.rs` writes them.

Verified byte-identical against the upstream v0.1.3 binary:

| input | output | result |
| --- | --- | --- |
| SIRV parts + 4 reads | 18 lines | **identical** |
| 144 221 Drosophila parts + 20 000 reads | **6 542 313 lines** | **identical** |

and the SIRV digest `3a2cdf42…` is the same one already recorded in the `aln-keep-temp` golden's
`seeds.txt.gz`, so it agrees with the reference too, not just with the binary.

The resulting uLTRA binary links only `libc++`, `libz` and `libSystem`, and produces correct NAMs with
**no namfinder anywhere on `PATH`**.

What this buys, beyond not writing a strobemer implementation:

- **It fixes Finding 2 at the root.** namfinder missing from `osx-arm64` is the only reason the
  bioconda recipe cannot be solved there. The port has no such dependency.
- **It halves the external-tool surface.** Only `minimap2` remains, and deferred improvement D2 is
  about removing that too.
- **Finding 8 dissolves.** There is no `os.system` shell string left to quote badly.

Cost: a C++17 compiler and zlib at **build** time, nothing at run time. `libz` is currently a dynamic
link; making it static is a Stage 7 distribution question, not a correctness one.

### Finding 27 — the reference's alignments are not reproducible across index builds

Two indexes built from the same GTF on the same machine, differing only in `PYTHONHASHSEED`, produce
**different alignments**. Measured on sirv-10k, 10 000 real ONT reads, with the *alignment* runs both
pinned to `PYTHONHASHSEED=0` so the only variable is the index:

| file | identical? |
| --- | --- |
| `refs_sequences.fa` | **no** |
| `seeds.txt.gz` | **no** (same content, different line order) |
| **`reads.sam`** | **no — 2 lines, 1 read of 9 996** |

One read in ten thousand, and it is not a subtle difference: that read flips between **minimap2's
alignment and uLTRA's**, 22 SAM fields against 14, `NM/ms/AS` against `XA/XC:NIC_novel/NM`. The seed
order changed uLTRA's alignment slightly, which tipped `output_final_alignments`' choice of which
aligner won.

There are **two** independent causes, and fixing one is not enough:

1. **`refs_sequences.fa` is written in dict-insertion order** (`uLTRA:341`), which for the flank half
   is Python **set**-iteration order. Different order → namfinder reports tied NAMs in a different
   order → `seeds.txt.gz` differs.
2. **`parts_to_segments` and `gene_to_small_segments` hold their triples in set-iteration order**,
   which survives into `segment_hit_locations` as Finding 21's correction explains.

Fixing only (1) leaves `reads.sam` differing; fixing both makes it identical. Measured:

| | `refs_sequences.fa` | `seeds.txt.gz` | `reads.sam` |
| --- | --- | --- | --- |
| unpatched | differs | differs | **differs** |
| sort `refs_sequences.fa` only | same | same | **still differs** |
| sort both | same | same | **identical** |

**Why a user would ever hit this.** They would not set `PYTHONHASHSEED` at all, so every index build
picks a fresh order. Rebuild your index — same GTF, same genome, same version — and roughly one read
in ten thousand aligns differently. On a 10 M read dataset that is ~1 000 reads, silently, with
nothing in the output indicating which run you are looking at.

**The fix is two small changes, both semantically free**: write `refs_sequences.fa` sorted by the
unpacked `(chr, start, stop)` key, and sort the two arrays at index build. The consumers de-duplicate,
so sorting removes nothing; it only removes the arbitrariness.

It **does** change results — 2 lines on this corpus against a `PYTHONHASHSEED=0` baseline — so it
lands on its own branch and the goldens are re-recorded against it, exactly as the NGSpeciesID port
had to do with `--seed`. Goldens recorded against the unfixed reference would pin one arbitrary
ordering out of many.

For the port this is load-bearing, not optional: the port sorts (its structures are `BTreeMap`), so
without this fix it could only ever match a reference index that happened to have been built with a
compatible hash seed.

### Finding 28 — there are only three live aligner call sites, and two library wrappers that nothing calls

The aligner surface is much smaller than the dependency list suggests. Auditing every function in
`help_functions.py` that touches parasail or edlib, against its callers:

| wrapper | engine | callers |
| --- | --- | --- |
| `help_functions.edlib_alignment` | `edlib.align(task="path")` | 4 |
| `help_functions.parasail_alignment` | `sg_trace_scan_16` → `_32` | 1 |
| `classify_read_with_mams.edlib_alignment` | `edlib.align(...)` | 3 |
| `help_functions.ssw_alignment` | `parasail.ssw` | **0 — dead** |
| `help_functions.parasail_local` | `parasail.sw_trace_scan_16` | **0 — dead** |

So the port needs `edlib.align` and one parasail entry point, and nothing else.

Measured call mix on sirv-10k, recorded from a real run:

| site | calls |
| --- | --- |
| `edlib.align`, all `mode="HW", task="path"` | 11 810 |
| `parasail_alignment` | 190 |

**The ratio is misleading, and worth stating so nobody optimises the wrong one.** edlib dominates the
*count* because it is used for every MAM-level accuracy check, but parasail is the DEFAULT branch of
`get_exact_alignment` — edlib only takes over above 20 kb — so parasail produces most of the final
CIGARs that reach `reads.sam`.

**Both are linked, not reimplemented**, for the reason Finding 26 gives: edlib is vendored as a single
1 482-line `.cpp` (MIT, 88 KB), and parasail comes from `libparasail-sys`, which bundles its C
sources. This settles the one thing the API does not define — which optimal location edlib returns
first, given that `help_functions.edlib_alignment` reads exactly `locations[0]`.

Verified by replaying calls recorded from the reference, comparing full outputs rather than scores:

| oracle | records | compared | result |
| --- | --- | --- | --- |
| `tests/edlib_oracle.rs` | **1 866** | edit distance, full locations list, CIGAR | identical |
| `tests/parasail_oracle.rs` | **194** | score, CIGAR | identical |

The edlib replay also answers a question that is not discoverable from Python: `python-edlib`
1.3.9.post1 does not expose which edlib it bundles, and the replay confirms it behaves as v1.2.7.

Recording these required `bench/recorder/sitecustomize.py` rather than ordinary monkey-patching:
uLTRA aligns inside worker processes spawned by `pc.Managers`, and under `spawn` those children
re-import everything, so a patched parent records nothing. `sitecustomize` is imported during every
interpreter's site initialisation, children included.

### Finding 29 — the n log n chaining path almost never runs, which is exactly why it needs an oracle

`align.py` picks between two chaining implementations per chromosome per read:

```python
if len(all_mems_to_chromosome) < 90:
    solutions, v = colinear_solver.read_coverage(all_mems_to_chromosome, max_allowed_intron)
else:
    solutions, v = colinear_solver.n_logn_read_coverage(all_mems_to_chromosome)
```

Measured by recording every call on real data:

| corpus | `read_coverage` | `n_logn_read_coverage` | largest mem set |
| --- | --- | --- | --- |
| SIRV, 10 000 ONT reads | 281 | **0** | 38 |
| Drosophila, 20 000 reads | 23 169 | **45** | 340 |

So the n log n branch is **0.19 %** of chaining calls on Drosophila and **never fires at all** on
every SIRV corpus — which is five of the seven registered corpora, including the one the whole
project developed against.

That is the argument for the oracle in one line. A branch that runs on two reads in a thousand is a
branch whose bugs reach production: it is too rare to notice in a spot check and common enough to
corrupt a large run. It also carries the port's most delicate code — a range-max segment tree whose
tie-breaks (`max(sorted(V, key=j_max, reverse=True), key=Cj)`, in three separate places) decide which
of several equal-scoring chains is returned.

`bench/oracle/chaining_calls.jsonl.gz` therefore keeps **all 45** recorded n log n calls and only a
sample of the far more numerous quadratic ones. `tests/chaining_oracle.rs` replays both and compares
the **full solution set** — every mem of every optimal chaining, in order — not the score, because a
score-only check passes a port that returns a different equally-scoring chain, which is precisely
what a mis-ordered tie-break produces.

Result: **445 recorded calls replay identically, 45 of them through the n log n path**, up to a
340-mem instance.

Corpus consequence: droso-20k is the only registered corpus that exercises this code at all. If it is
ever dropped to save time, the n log n path silently stops being tested.

### Finding 30 — the oracle harness was silently corrupting its own floats

Not a finding in the reference; a finding in **this port's verification**, and the most useful kind.

`read_coverage_mam_score` scores in floating point, so its oracle compares values **exactly**, with no
epsilon. 57 of 368 recorded calls failed by about `2e-13`, with the chosen chain identical every
time — the classic signature of reassociated arithmetic.

It was not the arithmetic. A verbatim Python re-implementation of the reference algorithm reproduced
all 368 recorded values exactly, which located the fault in the Rust side. Instrumenting both to dump
the `C` vector found the first divergence, and the input at that step read:

| | |
| --- | --- |
| Python | `v.val = 36.660000000000004` |
| Rust | `v.val = 36.66` |

from the *same* JSON file, whose text really does say `36.660000000000004`.

**`serde_json`'s default float parser is not correctly rounded.** Measured:

```
serde_json  36.660000000000004  ->  0x4042547ae147ae14   (== 36.66)
Rust literal 36.660000000000004 ->  0x4042547ae147ae15
```

One ULP, on every float in every oracle. The fix is the crate's `float_roundtrip` feature, which is
off by default; with it enabled all 368 calls replay identically under exact comparison.

**Why this is worth a numbered finding.** The obvious response to "differs by 2e-13" is to compare
with a tolerance. That would have passed, felt reasonable, and permanently blinded the oracle to real
1-ULP divergences — in a program where `argmax` over near-equal float scores decides which alignment
is reported. The harness would have been green and wrong, which is the failure mode this project
exists to avoid.

Two rules earned:

- **An oracle must round-trip its own data losslessly, and that must be tested.** Any float-carrying
  oracle added later inherits this; `float_roundtrip` is now on for the workspace.
- **When a comparison fails by an amount that looks like rounding, find the cause before loosening
  the comparison.** The epsilon is the last resort, not the first.

### Finding 31 — one chaining variant is unreachable on every corpus, and it contains a live bug

`classify_read_with_mams.main` picks between two MAM chainers on size:

```python
if len(mam_instance) > 200:
    mam_solution, value, unique = colinear_solver.n_logn_read_coverage_mams(mam_instance, overlap_threshold=5)
else:
    mam_solution, value, unique = colinear_solver.read_coverage_mam_score(mam_instance, overlap_threshold=20)
```

Measured by recording every call:

| corpus | `read_coverage_mam_score` | `n_logn_read_coverage_mams` | largest MAM set |
| --- | --- | --- | --- |
| SIRV, 10 000 reads | 368 | **0** | 100 |
| Drosophila, 20 000 reads | 21 365 | **0** | **171** |

**The >200 branch does not fire on any registered corpus**, and 171 is not close enough to 200 to be
luck — no real read in either dataset produces that many MAMs. Note also that the two branches use
*different* overlap thresholds, 5 against 20, so they are not interchangeable implementations of one
function; reaching the second changes the scoring.

Two consequences.

**It cannot be verified from real data, so its oracle is synthetic and labelled as such.** Real
recorded MAM sets are concatenated until they cross the threshold, re-sorted by `y` and re-indexed,
and the *reference* is run on them to record the expected output — 25 inputs of 217 to 421 MAMs. That
is weaker than a real corpus, because the concatenation may not produce MAM geometries a genome
actually yields, and `bench/oracle/mam_chaining_calls.jsonl.gz` marks those records `synthetic: true`
so the distinction cannot be lost.

**It contains a bug that has therefore never been hit.** The I-tree query returns a `j_max`, which is
**negative** for the padding leaves added to round the leaf count up to a power of two, and the code
then does:

```python
prev_end = max(mam.c - 1, mams[j_prime_b].d)
ovl_penalty = mams[j_prime_b].d - (mam.c - 1) + 0.0001 if mam.x != mams[j_prime_b].y else 0.0001
```

`mams[-1]` does not raise in Python; it is the **last** MAM. So on a padding hit the overlap penalty
is computed against an unrelated MAM at the other end of the list. The sibling
`n_logn_read_coverage` avoids this by only using `j_prime` as a traceback index, never to look up a
neighbour.

The port reproduces the wrap-around, because the alternative is a silent behaviour change in code
that is already hard to test. It is marked in `colinear.rs` as a faithful bug, and it is a candidate
for *Deferred improvements* once the port is exact.

> Findings 32+ will be added as the port proceeds. The NGSpeciesID port accumulated 30 and they were
> the most useful artifact of the project.

---

## Part 7 — Repository slimming (goal 6)

Measured on a full mirror clone.

| | |
| --- | --- |
| repository on GitHub | **2 143 708 KB ≈ 2.1 GB** |
| unique blob bytes across all history | **2 234.6 MB**, 1 306 unique blobs, 3 429 objects, 1 pack |
| **working tree at HEAD** | **7.5 MB**, 68 files |

Unique blob bytes by top-level directory — **deduplicated by object id**, so this is not the
double-counting trap that NGSpeciesID hit as its Finding 21:

| bytes | path |
| --- | --- |
| **2 217.5 MB** | **`data/`** |
| 7.2 MB | `evaluation/` |
| 5.4 MB | `modules/` |
| 2.2 MB | `uLTRA` |
| 1.3 MB | `torkel2/` |
| 0.3 MB | `README.md`, `test/` |
| 0.2 MB | `torkel/` |
| < 0.1 MB | `setup.py`, `scripts/`, `.travis.yml`, `INSTALL.sh`, `travis.yaml` |

**`data/` is 99.2 % of the repository.** And — unlike NGSpeciesID — **zero blobs live at more than one
path**, so `git rev-list --objects` builds a correct removal list directly and the tooling is simpler
here than it was there.

The largest single items, all historical and none present at HEAD:

| bytes | path |
| --- | --- |
| 264.7 MB | `data/results/pacbio_alzheimer_success_cases_sam_files.tar.bz2` |
| 163.1 MB | `data/results/sim.csv.a.bz2` |
| 147.7 MB | `data/results/sim.csv.b.bz2` |
| 140.0 MB | `data/all_transcripts_ENSEMBL.txt.bz2.part-aa` |
| 5 × 93.1 MB | `data/1M_NIC.fa.bz200` … `bz204` |

`data/` at HEAD is only 4.0 MB / 26 files, **no code references it**, and the README links to exactly
one part of it — `data/images`, for the small-exon examples. So the reviewed removal list is
"everything ever under `data/` except the images at HEAD", and the repository goes from 2.1 GB to
roughly 17 MB.

The rewrite is a force-push. It is reviewed on its own, on its own branch, and it is **not** performed
without an explicit go-ahead.

---

## Part 8 — The staged plan

The rule from the template: **build the equivalence harness before porting anything**, then port one
stage at a time and keep it green.

### Stage 0 — harness first, no Rust

- `bench/setup_reference_env.sh` — the one conda line above, pinned.
- `bench/corpora.tsv` + `bench/corpora_fetch.local.sh` — SIRV real (10k/100k/full), SIRV simulated
  (err0, err7%), Drosophila, GENCODE. Pinned sha256s; the reads themselves are **not committed**,
  because committing gigabytes into a repository whose stated goal is to shrink would be absurd.
- `bench/equivalence.sh` — carried across from `NGSpeciesID/bench/`, with its
  `record` / `verify` / `stage` / `cli` / `stable` / `cli_audit` subcommands, its
  `check_bin_fresh` staleness gate (rule 4) and its explicit `REF_PYTHON` / `PORT_BIN` resolution
  (rule 3). **`PYTHONHASHSEED=0` on every reference run**, per Part 5.
- `bench/dump_reference.py` — stage oracles. For uLTRA these are not optional: the index stage has no
  byte-identity contract at all, so the canonical text rendering of the 14 index structures **is**
  the contract. Also needed: recorded parasail calls, recorded edlib calls, recorded namfinder
  seed parses.
- Record goldens on ≥5 corpora. Prove the harness has teeth by running deliberately broken "ports"
  against it.

### Stage 1 — CLI contract

Three subcommands (`pipeline`, `index`, `align`), **22 unique long flags** plus `--version`/`-h` and
four positionals; `pipeline` takes 24 `add_argument` calls, `align` 19, `index` 9. The
`--ont`/`--isoseq` presets are applied **after** parsing and silently override an explicit `--s`
(Finding 12). Hand-written argparse-compatible parser, exit codes and stderr captured as goldens.
Everything past validation exits non-zero until implemented.

### Stage 2 — `index` — **DONE, and it was the biggest single win**

GTF parser replacing gffutils and its sqlite; the augmented-gene graph; the segment, exon, flank and
splice structures; sequence extraction. Verified against the stage oracle, not against bytes:
**119 of 120 structures match across all six corpora, with one approved divergence** (Finding 22) and
zero failures.

Measured, port against reference, same machine:

| corpus | reference (default) | reference (`--disable_infer`) | **port** |
| --- | --- | --- | --- |
| SIRV | 0.32 s, 62.9 MB | — | **0.08 s, 2.6 MB** |
| Drosophila | 240.3 s, 863.1 MB | 21.1 s, 827.0 MB | **1.12 s, 828.3 MB** |

**215× against the default path and 19× against `--disable_infer`**, and no 399 MB sqlite database is
written at all. On SIRV the memory drop is 24×, which is mostly the Python interpreter and gffutils
not being there.

Two honest caveats. The Drosophila memory is essentially unchanged because both implementations hold
the genome and every extracted sequence in memory at once; streaming that is a separate change, not
something the rewrite gave for free. And the port currently renders the index rather than writing a
loadable on-disk format — index *construction* is what this stage verifies, and the format arrives
when `align` needs to read one.

Bugs the reference has here, all found by running the new corpus rather than by reading:
Findings 13, 14, 22, 23, 24, 25.

### Stage 3 — seeds and I/O — **seeding DONE; read I/O outstanding**

**Seeding is done and exact**, by vendoring and linking namfinder rather than reimplementing it
(*Finding 26*). Byte-identical to the upstream binary on 4 reads and on 6 542 313 output lines of
Drosophila, and equal to the digest already recorded in the `aln-keep-temp` golden.
`namfinder::argv_for` reproduces `find_nams_namfinder`'s argument derivation including the
`(strobe_size + 1)//3` and `//5` thinning arithmetic, which is truncating integer division and must
stay so. Finding 8 has dissolved: there is no `os.system` string left.

**Still outstanding in this stage:** uLTRA's own read handling — the `reads_tmp.fa.gz` writer, the
gzip seed reader, and the read parser with a real error on malformed fastq instead of the silent
`zip()` truncation of *Finding 6*.

### Stage 4 — the alignment core — **aligners DONE; chaining outstanding**

**The aligners are done and exact.** Only three live call sites exist (*Finding 28*); edlib is
vendored and parasail comes from `libparasail-sys`, and both replay recorded reference calls
identically — 1 866 and 194 records respectively, comparing CIGARs and locations rather than scores.

**The colinear solver is done and exact.** Both implementations — `read_coverage` (quadratic) and
`n_logn_read_coverage` with its range-max segment tree — replay recorded reference calls identically:
445 calls, 45 of them through the n log n path, compared on the full solution set rather than the
score. See *Finding 29* for why the rare branch got the careful treatment.

**The MAM layer's chaining is done**: `read_coverage_mam_score` replays 368 recorded calls
identically, including the float score compared **exactly** — see *Finding 30* for what that
uncovered.

**Both MAM chainers are ported.** `read_coverage_mam_score` against 368 real recorded calls and
`n_logn_read_coverage_mams` against 25 synthetic ones (*Finding 31*): 393 total, all identical.

**Still outstanding:** `add_segment_to_mam` and the two `get_unique_*` helpers,
`classify_read_with_mams.main`'s orchestration, `align.annotate_guaranteed_optimal_bound`, and
`align.find_exons`.

### Stage 5 — SAM output and the minimap2 merge

`sam_output.py` and `prefilter_genomic_reads.py`. `reads.sam` byte-identical against the goldens.

### Stage 6 — parallelism

One shared immutable index, threads not processes (Finding 4). This is where goal 4 is won; verify
that peak RSS is flat in `--t` rather than +337 MB per worker.

### Stage 7 — distribution

cargo-zigbuild to `x86_64-unknown-linux-gnu.2.17` and aarch64, plus a bioconda recipe. Under
decision B the recipe's only runtime dependency is `minimap2`, which is present on all four subdirs
— so the port's own install works on Apple Silicon without waiting on namfinder#1. Test the way
a stranger installs: a real `pip`/`conda` install into a clean environment, and a release binary
downloaded through a **browser** so macOS actually sets the quarantine xattr — never via
`gh release download` (rule 1).

### Running in parallel, off `master`, not waiting for the port

A `fix/installation` branch: Findings 1, 7 and 9, plus the `.gitignore` the repository has never had.
And **[ksahlin/namfinder#1](https://github.com/ksahlin/namfinder/issues/1)** — filed — for the
osx-arm64 bioconda build (Finding 2), which fixes the reference's install on Apple Silicon.

---

## Deferred improvements

These land **after** the port is exact, each in its own commit, each measured against the byte-identical
port as the control. That ordering is the whole point: none of them can be evaluated honestly until
there is a known-good baseline to compare against.

### D1 — multi-context seeds (MCS), as in strobealign

uLTRA seeds with randstrobes (two strobes), and a seed is found only if **both** strobes match.
strobealign v0.15.0 changed its index layout so that a miss can fall back to looking up a single
strobe — a *partial seed* — and enabled this by default in v0.16.0. Two strategies exist there: the
default searches all full seeds and falls back to partial only if *none* hit, while `--mcs` falls
back **per seed** and is more accurate and slower.

Reported to improve mapping rate and accuracy for reads up to ~200 nt
([Tolstoganov et al., Genome Biol 2026](https://doi.org/10.1186/s13059-026-04017-x)). uLTRA's reads
are long transcriptomic reads, i.e. far above that range, so **the gain here is a hypothesis to test,
not a result to import** — the published benefit is at a read length uLTRA does not operate at. What
may still transfer is sensitivity in short exons, which is exactly where uLTRA claims its advantage,
and where a full randstrobe may not fit. That is the experiment: accuracy per exon size, MCS on
versus off, against the exact port.

Prerequisite: decision B, since partial-seed lookup requires uLTRA to own the index layout rather
than parse namfinder's output.

### D2 — drop minimap2 by seeding against the full genome as well as the annotation

Today minimap2 is 33 % of default wall clock and one of two external binaries. It does **two** jobs,
and both must be replaced for it to go away:

1. **Genomic-read detection.** `prefilter_genomic_reads` maps with minimap2, builds an interval tree
   of uLTRA-indexed regions, and calls a read genomic when more than `--genomic_frac` (0.1) of its
   aligned length falls outside them. Replacing this means indexing the **whole genome**, not just
   annotated regions plus `--flank_size` flanks — a much larger index, which interacts directly with
   goal 4.
2. **Primary-alignment cross-check.** `output_final_alignments` keeps the better of minimap2's and
   uLTRA's alignment per read. On the 4-read smoke corpus minimap2 still wins occasionally — the tool
   itself reports *"1 primary alignments had slightly better score with alternative aligner (typically
   ends bonus giving better scoring in ends, which needs to be implemented in uLTRA)"*. So this is an
   **accuracy** question, not only an engineering one, and the ends-bonus gap is a concrete,
   already-documented prerequisite.

The payoff is large and should be stated plainly: it removes the last external dependency, making the
binary genuinely self-contained, **and it invalidates the 2.5× ceiling recorded in Part 2** — that
ceiling exists only because 33 % of wall clock belongs to a subprocess a port cannot speed up. Note
the replacement is not free: seeding against the full genome costs time and memory that minimap2
currently absorbs, so the ceiling moves rather than vanishes, and by how much is a measurement.

### D3 — MEM seeding instead of NAMs

Carried over from *The seeding decision*. MEMs may be more sensitive than NAMs; changing the seeds
changes every output file, so it needs the exact port as a control. Also requires reviving
`get_mem_records`, which is dead *and* broken (Finding 7).

### D4 — drop `dill` from the Python reference

Four `defaultdict(lambda: array("L"))` factories are the only reason dill is a dependency
(Hypothesis 3). Module-level named functions remove it. This is an upstream commit against the
Python tool, independent of the port, and it changes the index pickles — so it must not land between
recording goldens and verifying against them.

## Working agreements

- Python is the spec. Divergences are numbered Findings, never silent.
- Every number is read off a run. Nothing is derived by reimplementing the tool's logic.
- No tagging, releasing or publishing without the author asking for it.
- No dataset is committed or published unless the author has confirmed it is public.
- The history rewrite is not force-pushed without an explicit go-ahead.
