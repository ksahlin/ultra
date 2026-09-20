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
set. The default path infers genes and transcripts and is slower still.

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

**Decision taken: B.** Strobemers are computed inside the Rust binary; namfinder stops being a
runtime dependency. **A is the fallback** if the NAM oracle cannot be made exact — that is a
measurement, not a preference, and it is made before the seeding stage is written, not after.

C is *Deferred improvements*. MEMs may well be more sensitive than NAMs, but changing the seeds
changes every output file, so C cannot be evaluated until the port is exact enough to serve as the
control. Port first, then measure MEMs against a known-good baseline. Reviving it would also mean
reviving `get_mem_records`, which is dead *and* broken (Finding 7).

Because B is a byte-identity risk concentrated in one place, it gets the treatment spoa got in the
NGSpeciesID port: `bench/dump_reference.py --stage seeds` records namfinder's NAM output on every
corpus, and `rust/tests/nam_oracle.rs` replays it. The port does not proceed past seeding until that
oracle is green or A has been adopted with a Finding explaining why.

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

> Findings 13+ will be added as the port proceeds. The NGSpeciesID port accumulated 30 and they were
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

### Stage 2 — `index`, and it is the biggest single win

GTF parser replacing gffutils and its 3.7 GB sqlite; the augmented-gene graph
(`create_augmented_gene.py`, 641 lines); the segment/exon/flank structures. Verified against the
stage oracle, not against bytes. **Expected: 140 s and 3.7 GB of sqlite → seconds and no database.**

### Stage 3 — seeds and I/O

**Strobemers computed natively (decision B).** The NAM oracle comes first: record namfinder's output
across every corpus, then make the Rust seeder reproduce it exactly. `--s` and `--thinning` map onto
namfinder's `-k/-s/-l/-u` exactly as `find_nams_namfinder` derives them, including the
`(strobe_size + 1)//3` and `//5` thinning arithmetic, which is integer division and must stay so.

Also here: the read parser, with a real error on malformed fastq instead of silent truncation
(Finding 6). Finding 8 dissolves — there is no `os.system` line left to quote badly.

### Stage 4 — the alignment core

`colinear_solver.py` + `range_query_max_search_tree.py` + `classify_read_with_mams.py` +
`align.py` — the 59 % of wall clock. Each aligner call site gets a recorded oracle **before** it is
written, the way the NGSpeciesID port measured edlib's tie-break rather than guessing it.

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

## Working agreements

- Python is the spec. Divergences are numbered Findings, never silent.
- Every number is read off a run. Nothing is derived by reimplementing the tool's logic.
- No tagging, releasing or publishing without the author asking for it.
- No dataset is committed or published unless the author has confirmed it is public.
- The history rewrite is not force-pushed without an explicit go-ahead.
