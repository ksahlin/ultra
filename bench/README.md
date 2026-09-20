# The uLTRA equivalence harness

The Python implementation in `../uLTRA` and `../modules/` is the **reference**. This directory
exists to answer one question mechanically: *does the Rust port produce the same bytes?*

It was built **before** any porting started, which is the whole point. A harness written after the
port exists gets shaped, unconsciously, around what the port already does.

## Quick start

```sh
bench/setup_reference_env.sh          # one conda line; builds namfinder where bioconda lacks it
conda activate ultra_ref
bench/equivalence.sh check            # ALWAYS run this first
bench/equivalence.sh record           # run the reference, write goldens
bench/equivalence.sh stable           # prove the goldens are reproducible
bench/equivalence.sh verify           # run the port against them
```

## The contract, and where it is unusual

Most of the port is ordinary byte-identity work. Two things are not, and both are measured in
`../PORTING.md` Part 5 rather than assumed here.

**`reads.sam` is a byte-identity target.** It carries `@SQ` lines and nothing else — no `@PG`, no
`@HD`, so no version string, no command line, no absolute paths. This was checked, because it is
exactly the trap the NGSpeciesID port hit with medaka's BAM.

**The index is not, and never can be.** It is 20 Python pickles plus a gffutils sqlite database. A
Rust port will not produce those, so byte-identity is not a coherent goal for them. They are
classified `oracle` and covered by a semantic stage oracle instead.

`equivalence.sh classify()` is the whole contract in one function:

| class | applies to | why |
| --- | --- | --- |
| `hash` | `reads.sam`, `*.fastq`, everything else | byte-for-byte |
| `gunzip` | `*.gz` | gzip writes an **mtime** into header bytes 5–8, so `seeds.txt.gz` differs every run while its contents are identical. Hash the decompressed stream. |
| `scrub_pg` | `minimap2.sam`, `indexed.sam`, `unindexed.sam` | these carry minimap2's `@PG` with its version and the verbatim command line, including whatever paths were passed. Only visible under `--keep_temporary_files`. |
| `exists` | `*_errors.*`, `*_stderr*`, `*.stderr` | content is timings. Record presence and size; a missing log is still a failure. |
| `oracle` | `*.pickle`, `database.db` | see above |

## `PYTHONHASHSEED=0` is not optional

Two identical reference runs differ in **17 of 34** output files by default and **3 of 34** with
`PYTHONHASHSEED=0`. The 14 it fixes are index pickles and `refs_sequences.fa`, which serialise sets
and dicts whose iteration order is hash-seed dependent. The 3 that remain are the two gzip mtimes and
minimap2's timing log, all handled by `classify()`.

`reads.sam` is deterministic either way — but the harness pins the seed on every reference run
regardless, because a contract that holds by luck is not a contract.

## Corpora

`corpora.tsv` pins a sha256 for the reference, the annotation **and** the reads of every corpus.
`corpora_resolve.sh` verifies all three before a run and exits non-zero on a mismatch. No read data
is committed: this repository's other goal is to shrink from 2.1 GB.

```sh
bench/corpora_resolve.sh --list
bench/corpora_resolve.sh --check-all
bench/corpora_resolve.sh sirv-10k     # -> ref<TAB>gtf<TAB>reads
```

`$ULTRA_DATA` (default `$HOME/data`) locates the public corpora; URLs are fetched into
`$ULTRA_BENCH_CACHE` (default `$HOME/.cache/ultra-bench`).

**`droso-20k` is not optional.** Six of the seven corpora share one 7-sequence reference and a
339-line GTF, which leaves the index stage — the largest rewrite in the port — effectively
unexercised. Adding a real genome and a real annotation found `PORTING.md` Finding 13 within minutes:
the reference *crashes* on a standard Ensembl GTF.

## Has the harness got teeth?

Measured, by building deliberately-broken "ports" and checking the harness rejects them. On the smoke
corpus, 9 cases:

| fake port | what it breaks | result |
| --- | --- | --- |
| `p0_faithful` | nothing — the control | **9 ok, 0 failed** |
| `p1_drop_line` | deletes the last SAM line | 8 caught |
| `p2_flip_base` | flips one base in one `SEQ` | 7 caught |
| `p3_preset_bug` | ignores `--ont`/`--isoseq` (Finding 12) | 3 caught |
| `p4_empty` | exits 0, writes nothing | 9 caught |
| `p5_gz_content` | changes `.gz` content, not its mtime | 8 caught |

`p5` matters more than it looks: it confirms the `gunzip` rule compares **content** and has not
quietly degenerated into a no-op that passes everything.

### What the teeth test found out about the corpora

`p3_preset_bug` was caught by `aln-ont`, `aln-ont-s12` and `aln-isoseq-s12` — but **not** by
`aln-isoseq`. On four reads, `--isoseq` is indistinguishable from no preset at all, because it sets
`s = 10`, which is already the default, and the only other thing it changes (`min_acc` 0.5 → 0.8)
does not alter the output of a 4-read corpus.

That is `PORTING.md` Finding 10 caught in the act. Re-running the *same* fake port against
`sirv-10k` — 10 000 real ONT reads, same reference, same cases — catches all three:

```
  ok      idx-default
  FAIL    aln-ont
  FAIL    aln-isoseq          <- passed on smoke, fails here
  FAIL    aln-isoseq-s12
  ======== verify: 1 ok, 3 failed, 0 skipped, of 4 ========
```

So the argument for the `core` tier is not a matter of taste: **a green smoke run means nothing on
its own**, and here is a concrete bug it would have shipped.

## Layout

| file | role |
| --- | --- |
| `setup_reference_env.sh` | pins the reference environment; builds namfinder from source on subdirs bioconda does not cover |
| `corpora.tsv` | corpus registry, sha256 for every input |
| `corpora_resolve.sh` | resolves a corpus to verified local paths |
| `cases.tsv` | the case matrix — 34 invocations |
| `equivalence.sh` | `check` / `record` / `verify` / `stable` / `list` |
| `golden/<corpus>/<case>/` | `manifest.tsv` and `exit` |
