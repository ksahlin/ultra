#!/usr/bin/env python3
"""Generate rust/src/text/*.txt from the recorded CLI goldens.

The port's fixed strings are EXTRACTED, never retyped. Retyping a 700-byte
argparse help block by hand and expecting it to be byte-identical is not a
plan; and if the reference's help changes, re-running this is the diff.

Run from the repository root:  python3 bench/extract_cli_text.py
"""
import pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
G = ROOT / "bench" / "golden" / "cli"
OUT = ROOT / "rust" / "src" / "text"

def read(case, stream="stdout"):
    p = G / case / stream
    if not p.exists():
        sys.exit(f"missing golden: {p} -- run 'bench/equivalence.sh cli record' first")
    return p.read_text()

def usage_block(text):
    """The wrapped usage block.

    In --help output it is followed by a blank line. In error output there is
    no blank line -- the "<prog>: error: ..." line follows immediately -- so
    both terminators are honoured. (Getting this wrong is what the cross-check
    below caught the first time this script was run.)
    """
    out = []
    for line in text.splitlines():
        if not line.strip():
            break
        if line.startswith("uLTRA") and ": error:" in line:
            break
        out.append(line)
    return "\n".join(out) + "\n"

OUT.mkdir(parents=True, exist_ok=True)
written = []

# Full help bodies (stdout of --help).
for case, name in [("no-args", "help_top"),
                   ("sub-help-pipeline", "help_pipeline"),
                   ("sub-help-index", "help_index"),
                   ("sub-help-align", "help_align")]:
    body = read(case)
    (OUT / f"{name}.txt").write_text(body)
    written.append((f"{name}.txt", len(body)))

# Usage blocks, taken from the SAME goldens so they cannot drift apart.
for case, name in [("no-args", "usage_top"),
                   ("sub-help-pipeline", "usage_pipeline"),
                   ("sub-help-index", "usage_index"),
                   ("sub-help-align", "usage_align")]:
    u = usage_block(read(case))
    (OUT / f"{name}.txt").write_text(u)
    written.append((f"{name}.txt", len(u)))

# --version.
v = read("version")
(OUT / "version.txt").write_text(v)
written.append(("version.txt", len(v)))

# Cross-check: the usage argparse prints on an ERROR must equal the usage it
# prints in --help. If these ever disagree, one of them is being retyped.
checks = [("bad-t", "usage_align"), ("bad-subcommand", "usage_top")]
for case, name in checks:
    from_err = usage_block(read(case, "stderr"))
    from_help = (OUT / f"{name}.txt").read_text()
    if from_err != from_help:
        sys.exit(f"MISMATCH: usage in {case}/stderr differs from {name}.txt\n"
                 f"--- err ---\n{from_err}\n--- help ---\n{from_help}")

for n, b in written:
    print(f"  wrote {n:22s} {b:6d} bytes")
print(f"cross-checked {len(checks)} usage blocks against error output: consistent")
