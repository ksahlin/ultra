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

Full `uLTRA align` including the minimap2 and namfinder steps, on a 16-core machine (12 performance
+ 4 efficiency). minimap2 and namfinder are *the same code in both*, so they dilute every figure
here:

| corpus | mode | `--t` | Python | Rust | |
| --- | --- | --- | --- | --- | --- |
| sirv-10k | `--disable_mm2` | 8 | 4.3 s | **1.5 s** | 2.9x |
| sirv-10k | minimap2 | 8 | 6.0 s | **2.6 s** | 2.3x |
| sirv-10k | `--disable_mm2` | 3 | 7.2 s | **3.0 s** | 2.4x |
| sirv-10k | minimap2 | 3 | 10.6 s | **5.8 s** | 1.8x |
| droso-20k | `--disable_mm2` | 8 | 9.1 s | **6.7 s** | 1.4x |
| droso-20k | minimap2 | 8 | 14.7 s | **10.0 s** | 1.5x |
| droso-20k | `--disable_mm2` | 3 | 13.8 s | **9.3 s** | 1.5x |
| droso-20k | minimap2 | 3 | 20.5 s | **14.4 s** | 1.4x |

**The speedup depends on the annotation, so quote the range, not a single number**: 1.4x on
Drosophila, up to 2.9x on SIRV. Earlier versions of this file said "about 2x" on the strength of
SIRV alone, and at that point the port was in fact **2.3x slower** on Drosophila — see PORTING.md
*Finding 45*, which is also why SIRV could not show it.

**Peak memory does not grow with `--t`.** 1077 MB at one thread and 1078 MB at sixteen on sirv-10k,
where the Python implementation grows from 1139 MB to 2845 MB because it starts one process per
core and each re-loads the index. At 16 threads that is 2.6x less memory and 2.3x faster at the
same time.

Most of that peak is namfinder's seed index rather than uLTRA: on sirv-10k, namfinder alone accounts
for 1042 MB of the 1078 MB total.

How it was verified
-------------------

The Python implementation is the specification and byte-identical output is the acceptance
criterion. `reads.sam` is compared **read by read** against the Python implementation, and every
difference has to fall into one of the documented classes below — none is tolerated as noise.
Behind that sit oracle tests replaying recorded Python calls through the Rust code: 1866 edlib
alignments, 2182 prefilter decisions, 2081 CIGAR merges, 900 SAM records byte for byte, and others.

| corpus | mode | reference | port | identical | the **alignment** differs |
| --- | --- | --- | --- | --- | --- |
| sirv-100k | `--disable_mm2` | 100 000 | 100 000 | **99.88 %** | 36 (0.036 %) |
| sirv-100k | minimap2 | 99 934 | 100 000 | 92.08 % | 9 (0.009 %) |
| droso-20k | `--disable_mm2` | 20 000 | 20 000 | 86.04 % | 467 (2.34 %) |
| droso-20k | minimap2 | 18 310 | 20 000 | 92.45 % | 56 (0.31 %) |
| droso-200k | `--disable_mm2` | 200 000 | 200 000 | 86.16 % | 4 542 (2.27 %) |
| droso-200k | minimap2 | 176 828 | 200 000 | 90.93 % | 551 (0.31 %) |

The last column is the one that matters, and it is worth separating from the rest. Differences
split three ways:

| | droso-200k `--disable_mm2` | droso-200k minimap2 |
| --- | --- | --- |
| **A. different locus** (RNAME/POS) | 2 139 (1.07 %) | 441 (0.25 %) |
| **B. same locus, different CIGAR** | 2 104 (1.05 %) | 110 (0.06 %) |
| **C. same alignment, reporting only** | 20 143 (10.07 %) | 15 481 (8.76 %) |
| mapped on one side, unmapped on the other | 299 (0.15 %) | — |

**About 90 % of all differences are class C — the alignment is identical and only how it is written
down differs** (QUAL orientation, SEQ strand, `XA` ordering). Where the alignment genuinely differs,
identity against the genome is the same on both sides, median 97.8 %; the port simply aligns
slightly fewer bases, **-0.225 % of aligned bases overall**, concentrated in 2 % of reads.

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
