#!/usr/bin/env python3
"""Record every aligner call the reference makes, for replay against the port.

WHY. PORTING.md's plan for stage 4 says each aligner call site gets a recorded
oracle BEFORE it is written -- the way the NGSpeciesID port measured edlib's
tie-break rather than guessing it. There are exactly three live sites:

  1. classify_read_with_mams.edlib_alignment -> edlib.align(mode=HW, task=locations, k=...)
     the MAM accuracy check
  2. help_functions.edlib_alignment          -> edlib.align(task=path, mode=HW)
     the long-sequence alignment path
  3. help_functions.parasail_alignment       -> sg_trace_scan_16, falling back to _32
     the default alignment path

(`ssw_alignment` and `parasail_local` exist and are called from nowhere.)

Each record captures the inputs and the FULL output, so a replay can check not
just the score but the CIGAR and the chosen location -- which is where
tie-breaks hide.

Usage:
  bench/record_aligner_calls.py --ref R.fa --reads X.fq --index IDX --out calls.jsonl
                               [--limit N] [--] [extra uLTRA align args]
"""
import argparse, json, os, sys

ap = argparse.ArgumentParser()
ap.add_argument("--ref", required=True)
ap.add_argument("--reads", required=True)
ap.add_argument("--index", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--limit", type=int, default=0, help="stop recording after N calls (0 = all)")
ap.add_argument("rest", nargs="*")
args = ap.parse_args()

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import edlib
from modules import help_functions, classify_read_with_mams

out = open(args.out, "w")
state = {"n": 0}

def rec(kind, payload):
    if args.limit and state["n"] >= args.limit:
        return
    state["n"] += 1
    out.write(json.dumps({"kind": kind, **payload}, sort_keys=True) + "\n")

# --- site 1 + 2: the raw edlib.align, wrapped once so both sites are captured
_edlib_align = edlib.align
def edlib_align_rec(query, target, **kw):
    r = _edlib_align(query, target, **kw)
    rec("edlib.align", {
        "query": query, "target": target,
        "mode": kw.get("mode"), "task": kw.get("task"), "k": kw.get("k", -1),
        "editDistance": r.get("editDistance"),
        "locations": r.get("locations"),
        "cigar": r.get("cigar"),
        "alphabetLength": r.get("alphabetLength"),
    })
    return r
edlib.align = edlib_align_rec
help_functions.edlib.align = edlib_align_rec
classify_read_with_mams.edlib.align = edlib_align_rec

# --- site 3: parasail, recorded at the wrapper so the _16 -> _32 fallback and
#     the derived alignment strings are captured together
_parasail_alignment = help_functions.parasail_alignment
def parasail_rec(s1, s2, **kw):
    r = _parasail_alignment(s1, s2, **kw)
    read_aln, ref_aln, cigar_string, cigar_tuples, score = r
    rec("parasail_alignment", {
        "s1": s1, "s2": s2, "kw": {k: v for k, v in kw.items()},
        "read_aln": read_aln, "ref_aln": ref_aln,
        "cigar": cigar_string, "score": score,
    })
    return r
help_functions.parasail_alignment = parasail_rec
import modules.align as align_mod
align_mod.help_functions.parasail_alignment = parasail_rec

# --- run the reference's align step
sys.argv = ["uLTRA", "align", args.ref, args.reads, os.environ.get("ULTRA_RECORD_OUT", "/tmp/_rec_out"),
            "--index", args.index] + args.rest
os.makedirs(sys.argv[4] if False else os.environ.get("ULTRA_RECORD_OUT", "/tmp/_rec_out"), exist_ok=True)

ns = {"__name__": "__main__", "__file__": os.path.join(ROOT, "uLTRA")}
try:
    with open(os.path.join(ROOT, "uLTRA")) as fh:
        code = fh.read()
    exec(compile(code, os.path.join(ROOT, "uLTRA"), "exec"), ns)
except SystemExit:
    pass
finally:
    out.close()

n = sum(1 for _ in open(args.out))
print(f"recorded {n} aligner calls to {args.out}", file=sys.stderr)
