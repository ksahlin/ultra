The Rust port
=============

uLTRA was reimplemented in Rust. The command line is unchanged, so existing pipelines do not need
editing. This file is what a user needs: what got faster, what got smaller, and the ways the output
differs. The full engineering record — 42 numbered findings, every measurement, every deliberate
divergence and why it was taken — is in [PORTING.md](PORTING.md).

Installation
------------

| | Python | Rust |
| --- | --- | --- |
| runtime dependencies | parasail-python, pysam, dill, gffutils, intervaltree, edlib, namfinder, minimap2 | **minimap2** |
| Apple Silicon | **conda install impossible** — namfinder has no osx-arm64 build, so it must be compiled by hand | yes |
| artefact | a Python package plus an external namfinder binary | **one binary**, 650–870 KB |

namfinder, edlib and zlib are compiled into the binary. Prebuilt binaries for linux-x86_64,
linux-aarch64, macOS arm64 and macOS x86_64 are produced by `packaging/build-release.sh`; the Linux
ones need only glibc 2.17 (RHEL/CentOS 7, 2012).

Speed and memory
----------------

sirv-10k, full `uLTRA align` including the minimap2 and namfinder steps, on a 16-core machine
(12 performance + 4 efficiency):

| `--t` | Python | Rust | Python peak | Rust peak |
| --- | --- | --- | --- | --- |
| 3 *(default)* | 11.2 s | **6.2 s** | 1139 MB | **1077 MB** |
| 16 | 4.9 s | **2.2 s** | **2845 MB** | **1078 MB** |

Two things to read off it.

The speedup is about **2x end to end**, which understates the port, because minimap2 and namfinder
are *the same code in both* and account for roughly half the Rust version's runtime. Dropping the
minimap2 dependency entirely is a possible future change; it is what would move the rest.

**Peak memory does not grow with `--t`** — 1077 MB at one thread and 1078 MB at sixteen. The Python
implementation starts one process per core and each re-loads the index, so it grows to 2845 MB at
16 threads. At that thread count the Rust version uses **2.6x less memory and runs 2.2x faster at
the same time**.

On Drosophila (20 000 reads against BDGP6.46), the Rust version takes 100.7 s on one thread and
**22.4 s on eight**, memory flat at ~3.2 GB. Most of that peak is namfinder's seed index, not
uLTRA: on sirv-10k namfinder alone accounts for 1042 MB of the 1078 MB total.

How it was verified
-------------------

The Python implementation is the specification and byte-identical output is the acceptance
criterion. `reads.sam` is compared **read by read** against the Python implementation, and every
difference has to fall into one of the documented classes below — none is tolerated as noise.
Behind that sit oracle tests replaying recorded Python calls through the Rust code: 1866 edlib
alignments, 2182 prefilter decisions, 2081 CIGAR merges, 900 SAM records byte for byte, and others.

| corpus | mode | reference | port | identical | differ |
| --- | --- | --- | --- | --- | --- |
| sirv-100k | `--disable_mm2` | 100 000 | 100 000 | **99 879 (99.88 %)** | 121 |
| sirv-100k | with minimap2 | 99 934 | 100 000 | 92 020 (92.08 %) | 7 914 |
| droso-20k | `--disable_mm2` | 20 000 | 20 000 | 17 208 (86.04 %) | 2 792 |
| droso-20k | with minimap2 | 18 310 | 20 000 | 16 927 (92.45 %) | 1 383 |

**Every single difference in all four runs is accounted for by a documented cause** — checked by
decomposing each differing record into causes rather than bucketing it, so that a record differing
for two known reasons at once is not miscounted as unexplained.

Two things the table shows that a smaller corpus does not. The minimap2 runs are dominated by the
QUAL-orientation fix, which is why they look *less* identical than the `--disable_mm2` ones: on
sirv-100k, 7 905 of the 7 914 differences are that one fix. And the port emits more reads than the
reference in both minimap2 runs — 66 and 1 690 — which is the dropped-read defect below.

How the output differs
----------------------

**Reads that minimap2 cannot place are no longer dropped.** This is the one that can change a
result. `filter_reads_to_align` writes such a read to neither `indexed.sam` nor `unindexed.sam`, so
the merge that builds the final file has no slot for it, and `uLTRA:309` then moves that merged file
*over* uLTRA's own SAM — so the alignments are computed, written, and deleted. Measured: **1 690 of
20 000 (8.5 %) on Drosophila**, of which **1 283 had a real uLTRA alignment** and 407 were
unalignable; and 66 of 100 000 on SIRV, 57 of them aligned. The port keeps uLTRA's alignment where
there is one and uLTRA's `FLAG 4` record where there is not, so **every input read appears in the
output exactly once**. The branch handling this case exists in the Python source and can never run
(*Finding 39*).

**Mapped records carry the strand their FLAG says they do.** `align.py:550` writes every alignment
using the enclosing loop's `read_seq`, so when the chosen alignment is not the last orientation
examined, SEQ is the reverse complement of what its own FLAG and CIGAR describe. Checked against
the genome: the reference's SEQ scores **26.0 % median identity** to the locus it claims — what
unrelated sequence scores — against the port's **97.1 %**, and the port is the better match in
**181 of 181** records checked. 167 of 20 000 on Drosophila without minimap2 (*Finding 43*). This
matters more than it sounds: a `FLAG 4` record warns you not to trust its orientation, a mapped one
does not, so anything downstream reading sequence rather than re-deriving it gets wrong bases
silently.

**QUAL matches SEQ for reverse-strand reads.** The Python implementation restores SEQ to the read's
orientation but leaves QUAL in the aligned one. Only on the minimap2 path; 7 905 of 100 000 on SIRV
(*Finding 34*).

**Unaligned records have a stable sequence** — always the read as given, rather than whichever
orientation the strand search examined last (*Finding 38*).

**`XA:Z:` transcript lists are sorted.** The Python implementation joins a `set` of transcript ids,
so its order depends on `PYTHONHASHSEED` and is not reproducible between runs of the Python tool
itself (*Finding 33*).

**A few alignments differ where the Python implementation's choice was arbitrary**, because segments
tied on start coordinate come out in CPython's `set` iteration order. The port uses a total order
(*Finding 33*).

**The cost of that last one, stated plainly.** On Drosophila without minimap2, **21 reads of 20 000
(0.105 %) that the reference aligns are reported unaligned by the port**, and 2 go the other way.
It is not a seeding failure — at `--alignment_threshold 0` the port aligns 20 of the 21 — but a
marginally lower coverage crossing a threshold. SIRV shows zero in both directions at 10 000 and
100 000 reads, so the small corpus cannot see it (*Finding 44*).

Reproducing the comparison
--------------------------

`bench/equivalence.sh` runs the Python implementation and the Rust binary over the same corpora and
compares every output file. `bench/stage_diffs.tsv` lists the approved divergences — anything else
differing is a failure, and a listed divergence that has silently disappeared is also reported, so
the list cannot rot into a mute list.

```bash
bench/equivalence.sh check             # dependencies, corpora, reference
bench/equivalence.sh record sirv-10k   # run the Python implementation, store hashes
bench/equivalence.sh verify sirv-10k   # run the Rust binary, compare
```
