# Vendored namfinder

Source: https://github.com/ksahlin/namfinder
Tag:    v0.1.3
Commit: 7468624191e779a828ae549fdfcf684413df2b22
Licence: MIT (see LICENSE) -- same author as uLTRA.

Only `src/` and `ext/` are vendored. The upstream CMake build is not used;
`rust/build.rs` compiles the same source list that upstream's `salib` target
does, plus `main.cpp` with its `main` renamed, and `shim.cpp` which exposes
`run_strobealign` under a C ABI.

## Why vendored rather than shelled out to

namfinder is the only reason the bioconda recipe cannot be solved on
`osx-arm64` (PORTING.md Finding 2), and it is one of only two external
binaries uLTRA needs. Linking it removes both problems at once and needs
nothing on PATH at run time.

## Why linked rather than reimplemented

PORTING.md's seeding decision chose "compute strobemers inside the Rust
binary". Linking satisfies that and is **exact by construction**: the port
calls `run_strobealign`, the very function namfinder's own `main` calls, so
argument parsing, index construction, NAM finding and output formatting are
all namfinder's. Verified byte-identical against the v0.1.3 binary on 4 reads
and on 20 000 Drosophila reads against 144 221 parts (6 542 313 output lines).

A from-scratch randstrobe implementation would have had to reproduce all of
that from an oracle, for no benefit.

## Updating

1. Copy `src/` and `ext/` from the new upstream tag; update this file.
2. `cargo build --release`
3. `bench/equivalence.sh stage seeds verify` -- it compares against seeds
   recorded from the reference, so a behaviour change upstream shows up as a
   diff rather than as silence.
