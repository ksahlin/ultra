# Vendored edlib

Source:  https://github.com/Martinsos/edlib
Tag:     v1.2.7
Licence: MIT (see LICENSE)
Files:   edlib.cpp (1482 lines), edlib.h (277 lines)

## Why vendored

uLTRA makes 11 810 of its 12 000 aligner calls through edlib (measured on
sirv-10k), all of them `mode=HW, task=path`. The NGSpeciesID port found
`edlib_rs` and `rsedlib` both unusable and recorded vendoring the single .cpp
as the fallback; here it is also the *first* choice, for the same reason
namfinder is linked (PORTING.md Finding 26): exact by construction beats exact
by oracle.

What is NOT uniquely defined in edlib's API is which location it returns first
when several achieve the optimum, and `help_functions.edlib_alignment` reads
exactly `result['locations'][0]`. Linking the same implementation removes that
question entirely.

## Version

The reference environment uses `python-edlib` 1.3.9.post1, whose bundled C++ is
not labelled from Python. The version here was confirmed by REPLAY rather than
by reading: `rust/tests/edlib_oracle.rs` replays calls recorded from the
reference and compares edit distance, locations and CIGAR.
