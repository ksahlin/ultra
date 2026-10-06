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

The Python implementation is the specification, and byte-identical output was the acceptance
criterion. `reads.sam` is compared read by read against the Python implementation on real corpora,
and every difference is accounted for by one of the documented divergences below — not tolerated as
noise. Supporting this there are oracle tests that replay recorded Python function calls through
the Rust code: 1866 edlib alignments, 2182 prefilter decisions, 2081 CIGAR merges, 900 SAM records
compared byte for byte, and others.

On 10 000 real SIRV reads with `--disable_mm2`, **9990 of 10 000 records are byte-identical**, and
all 10 that differ fall into the classes below.

How the output differs
----------------------

**Reads that minimap2 cannot place are no longer dropped.** This is the one that can change a
result. `filter_reads_to_align` writes such a read to neither `indexed.sam` nor `unindexed.sam`, so
the merge that builds the final file has no slot for it and whatever uLTRA aligned is discarded.
The read vanishes. Measured: **4 of 10 000** on SIRV, and **1690 of 20 000 — 8.5% — on
Drosophila**, every one of them a read minimap2 could not place, which is the population uLTRA is
most likely to help with. The Rust version keeps uLTRA's alignment and prints the count. The branch
handling this case exists in the Python source and can never run; see PORTING.md *Finding 39*.

**Unaligned records have a stable sequence.** For a read nothing aligned, the Python implementation
writes whichever orientation the strand search happened to examine last, so SEQ is sometimes the
read and sometimes its reverse complement, with the FLAG giving no indication which. The Rust
version always writes the read as it was given. 2 of 10 000 on SIRV (*Finding 38*).

**QUAL matches SEQ for reverse-strand reads.** The Python implementation restores SEQ to the read's
orientation but leaves QUAL in the aligned one, so the two are reversed relative to each other.
This only arises on the minimap2 path. 759 of 10 000 records on SIRV (*Finding 34*).

**`XA:Z:` transcript lists are sorted.** When a read's splices match several transcripts, the
Python implementation joins a `set` of transcript ids, so the order depends on `PYTHONHASHSEED` and
is not reproducible between runs of the Python tool itself. The Rust version sorts. 4 of 10 000
(*Finding 33*, fifth site).

**A few alignments differ where the Python implementation's choice was arbitrary.** Segments tied
on start coordinate come out in CPython's `set` iteration order, which decides which is aligned
first and reaches the CIGAR. The Rust version uses a total order. 4 of 10 000 on SIRV, 0 of 2000 on
Drosophila. Reproducing it would mean reimplementing CPython's `setobject.c` probe sequence
(*Finding 33*).

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
