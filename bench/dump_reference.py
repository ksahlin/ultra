#!/usr/bin/env python3
"""Stage oracles: render the reference's internal state as canonical text.

WHY THIS EXISTS, and why it is not optional here.

For most stages of most ports, byte-identity of the output files IS the
contract. uLTRA's `index` stage cannot work that way: its output is 20 Python
pickles plus a gffutils sqlite database, and a Rust port will never produce
either. Byte-identity is not merely hard, it is incoherent.

So the contract for `index` is *semantic*: this script loads the reference's
pickles and emits a canonical, sorted, plain-text rendering of every structure;
the port emits the same rendering; the two are diffed. The rendering is the
specification.

Canonicalisation rules, chosen so that nothing machine- or run-dependent
survives:
  * dict / defaultdict  -> sorted by rendered key
  * set                 -> sorted
  * array('L')          -> flat [chr, start, stop] triples, SORTED. The stored
                           order is hash-seed dependent noise, not data - see
                           the note on --array-order below and PORTING.md
                           Finding 21.
  * bytes keys          -> decoded as three native unsigned longs and printed
                           as "chr/start-stop"

Usage:
    bench/dump_reference.py --index <indexdir> --out <dir>
    bench/dump_reference.py --index <indexdir> --out <dir> --only segment_to_gene
"""
import argparse
import os
import struct
import sys
from array import array

# array('L') is 8 bytes on LP64 (Linux, macOS) and 4 on Windows/32-bit, and the
# reference unpacks with the native 'LLL'. The port writes u64 little-endian,
# which is byte-identical on every platform uLTRA actually runs on. See
# PORTING.md Finding 20.
ITEMSIZE = array("L").itemsize
KEYLEN = 3 * ITEMSIZE


def decode_key(b):
    """A 3-long key -> 'chr/start-stop'."""
    if not isinstance(b, (bytes, bytearray)) or len(b) != KEYLEN:
        return repr(b)
    c, s, e = struct.unpack("LLL", bytes(b))
    return f"{c}/{s}-{e}"


ARRAY_ORDER = "sorted"   # set from --array-order


def render_value(v):
    if isinstance(v, array):
        # Flat [chr, start, stop] triples.
        #
        # ORDER. Only two structures are arrays -- parts_to_segments and
        # gene_to_small_segments -- and their stored order varies with
        # PYTHONHASHSEED, because it is the order of a set iteration upstream.
        # It is noise rather than data, established two ways:
        #
        #   * both consumers de-duplicate. gene_to_small_segments feeds a
        #     set(...) at classify_read_with_mams.py:266;
        #     parts_to_segments feeds segment_hit_locations, which is
        #     list(set(...)) then sorted at :244-245.
        #   * measured: reads.sam is byte-identical across hash seeds, while
        #     14 of the 20 pickles are not.
        #
        # So the port is free to build these in any order, and the oracle
        # sorts rather than pinning an order the reference does not really
        # have. --array-order stored restores the raw order for debugging.
        xs = list(v)
        trip = [(xs[i], xs[i + 1], xs[i + 2]) for i in range(0, len(xs) - 2, 3)]
        tail = xs[len(trip) * 3:]
        if ARRAY_ORDER == "sorted":
            trip.sort()
        out = ",".join(f"{c}/{s}-{e}" for c, s, e in trip)
        if tail:
            out += "|TRAILING:" + ",".join(map(str, tail))
        return out
    if isinstance(v, set):
        return ",".join(sorted(map(_atom, v)))
    if isinstance(v, dict):
        return ";".join(f"{_atom(k)}={render_value(x)}" for k, x in sorted(v.items(), key=lambda kv: _atom(kv[0])))
    if isinstance(v, (bytes, bytearray)):
        return decode_key(v)
    if isinstance(v, tuple):
        return "/".join(map(_atom, v))
    return _atom(v)


def _atom(x):
    if isinstance(x, (bytes, bytearray)):
        return decode_key(x)
    if isinstance(x, tuple):
        return "(" + ",".join(map(_atom, x)) + ")"
    if isinstance(x, float):
        return repr(x)
    return str(x)


def dump_mapping(o):
    """dict-like -> sorted 'key<TAB>value' lines."""
    items = [(_atom(k), render_value(v)) for k, v in o.items()]
    items.sort(key=lambda kv: kv[0])
    return [f"{k}\t{v}" for k, v in items]


def dump_set(o):
    return sorted(_atom(x) for x in o)


STRUCTURES = [
    "chr_to_id", "id_to_chr", "refs_lengths", "refs_id_lengths",
    "max_intron_chr",
    "segment_id_to_choordinates", "segment_to_ref", "segment_to_gene",
    "parts_to_segments", "gene_to_small_segments",
    "exon_choordinates_to_id", "flank_choordinates",
    "all_splice_pairs_annotations", "all_splice_sites_annotations",
    "splices_to_transcripts", "transcripts_to_splices",
    "ref_part_sequences", "ref_segment_sequences",
    "ref_exon_sequences", "ref_flank_sequences",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", required=True, help="an index folder built by the reference")
    ap.add_argument("--out", required=True, help="directory to write the rendering into")
    ap.add_argument("--only", action="append", default=None, help="restrict to these structures")
    ap.add_argument("--array-order", choices=("sorted", "stored"), default="sorted",
                    help="sorted (default) treats array order as noise; stored shows it raw")
    args = ap.parse_args()

    global ARRAY_ORDER
    ARRAY_ORDER = args.array_order

    # The index pickles embed a dill-serialised `defaultdict(lambda: array("L"))`
    # factory whose closure resolves through the `modules` package, so the
    # repository root must be importable or the load fails with
    # ModuleNotFoundError. Do it here rather than requiring a particular cwd.
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)

    import dill

    os.makedirs(args.out, exist_ok=True)
    wanted = args.only or STRUCTURES
    missing, written = [], []

    for name in wanted:
        path = os.path.join(args.index, f"{name}.pickle")
        if not os.path.exists(path):
            missing.append(name)
            continue
        with open(path, "rb") as fh:
            o = dill.load(fh)

        if isinstance(o, set):
            lines = dump_set(o)
        elif hasattr(o, "items"):
            lines = dump_mapping(o)
        else:
            lines = [render_value(o)]

        dest = os.path.join(args.out, f"{name}.txt")
        with open(dest, "w") as fh:
            fh.write("\n".join(lines))
            if lines:
                fh.write("\n")
        written.append((name, len(lines), os.path.getsize(dest)))

    for name, n, size in written:
        print(f"  {name:34s} {n:8d} lines  {size:10d} bytes")
    if missing:
        print(f"MISSING from {args.index}: {', '.join(missing)}", file=sys.stderr)
        return 1
    print(f"{len(written)} structures rendered into {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
