#!/usr/bin/env python3
"""Compare two reads.sam files read by read, and say WHAT differs.

    bench/compare_sam.py <reference.sam> <port.sam> <label> [genome.fa]

Splits differences into three tiers, because they mean very different things:

  A  different locus              RNAME or POS differs -- a different alignment
  B  same locus, different CIGAR  same place, different alignment detail
  C  same alignment, reporting    SEQ orientation, QUAL, tags only

plus mapped/unmapped disagreements and reads present on only one side.

A record can differ for two reasons at once -- a reversed QUAL and a reordered
XA, say. Bucketing by which FIELDS differ invents "unexplained" cases that are
really two known ones overlapping, so class C is decomposed by cause.

Given a genome, it also scores both sides' alignments under uLTRA's own scheme
(match +2, mismatch -2, gap open 3, extend 1, introns free) so a divergence can
be called better or worse rather than just different.
"""
import sys, collections

def load(p):
    d = collections.defaultdict(list)
    for line in open(p):
        if line.startswith('@'): continue
        f = line.rstrip('\n').split('\t')
        d[f[0]].append(f)
    return d

def tags(f):
    return {t[:2]: t[5:] for t in f[11:] if len(t) > 5}

ref, port, label = sys.argv[1], sys.argv[2], sys.argv[3]
r, p = load(ref), load(port)
only_p, only_r = sorted(set(p) - set(r)), sorted(set(r) - set(p))
common = [k for k in r if k in p]
diff = [k for k in common if r[k] != p[k]]
ident = len(common) - len(diff)

print(f"=== {label} ===")
print(f"  reference {len(r)} reads | port {len(p)} reads")
if only_r:
    print(f"  !! {len(only_r)} reads ONLY in the reference -- a BUG: {only_r[:4]}")
if only_p:
    al = sum(1 for k in only_p if not int(p[k][0][1]) & 4)
    print(f"  reads only in the PORT: {len(only_p)}  ({al} with an alignment, "
          f"{len(only_p)-al} FLAG 4)   [F39: the reference drops these]")
print(f"  shared {len(common)}: {ident} identical ({100*ident/max(len(common),1):.3f}%), {len(diff)} differ")
if not diff:
    print(); sys.exit()

tier = collections.Counter()
sub = collections.defaultdict(collections.Counter)
for k in diff:
    a, b = r[k][0], p[k][0]
    ua, ub = int(a[1]) & 4, int(b[1]) & 4
    if bool(ua) != bool(ub):
        tier['mapped on one side, unmapped on the other'] += 1
        sub['mapped on one side, unmapped on the other'][
            'reference aligns, port does not' if ub else 'port aligns, reference does not'] += 1
        continue
    if ua and ub:
        tier['C  both unmapped, record text differs'] += 1
        sub['C  both unmapped, record text differs']['SEQ/QUAL orientation [F38]'] += 1
        continue
    if a[2] != b[2] or a[3] != b[3]:
        tier['A  DIFFERENT LOCUS (RNAME/POS)'] += 1
        sub['A  DIFFERENT LOCUS (RNAME/POS)'][
            'different contig' if a[2] != b[2] else 'same contig, different POS'] += 1
    elif a[5] != b[5]:
        tier['B  same locus, different CIGAR'] += 1
        sub['B  same locus, different CIGAR']['segment tie order [F33]'] += 1
    else:
        tier['C  same alignment, reporting differs'] += 1
        c = sub['C  same alignment, reporting differs']
        ta, tb = tags(a), tags(b)
        if a[9] != b[9]: c['SEQ is the wrong strand in the reference [F43]'] += 1
        if a[10] != b[10] and a[9] == b[9]: c['QUAL orientation [F34]'] += 1
        if ta.get('XA') != tb.get('XA'): c['XA transcript order [F33]'] += 1
        if ta.get('XC') != tb.get('XC'): c['XC classification'] += 1
        if ta.get('NM') != tb.get('NM'): c['NM edit distance'] += 1
        if a[4] != b[4]: c['MAPQ'] += 1
        if a[1] != b[1]: c['FLAG'] += 1

order = ['A  DIFFERENT LOCUS (RNAME/POS)', 'B  same locus, different CIGAR',
         'C  same alignment, reporting differs', 'C  both unmapped, record text differs',
         'mapped on one side, unmapped on the other']
for t in order:
    if not tier[t]: continue
    pct = 100 * tier[t] / len(common)
    print(f"    {tier[t]:7d}  ({pct:6.3f}%)  {t}")
    for c, n in sub[t].most_common():
        print(f"            {n:7d}  {c}")
aln = tier['A  DIFFERENT LOCUS (RNAME/POS)'] + tier['B  same locus, different CIGAR'] \
      + tier['mapped on one side, unmapped on the other']
print(f"  --> the ALIGNMENT itself differs for {aln} of {len(common)} reads "
      f"({100*aln/len(common):.3f}%); the rest is reporting only")
print()

# ---------------------------------------------------------------------------
# Optional: score both sides against the genome, under uLTRA's own scheme.
# ---------------------------------------------------------------------------
if len(sys.argv) > 4:
    import re, statistics
    fa = sys.argv[4]
    cand = [k for k in r if k in p and not int(r[k][0][1]) & 4 and not int(p[k][0][1]) & 4
            and (r[k][0][2] != p[k][0][2] or r[k][0][3] != p[k][0][3] or r[k][0][5] != p[k][0][5])]
    if cand:
        need = {r[k][0][2] for k in cand} | {p[k][0][2] for k in cand}
        genome = {}; name = None; buf = []
        for line in open(fa):
            if line[0] == '>':
                if name in need: genome[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            elif name in need: buf.append(line.strip())
        if name in need: genome[name] = ''.join(buf)

        OPEN, EXT, MATCH, MIS = 3, 1, 2, -2
        def score(f):
            seq, cig, ctg, pos = f[9], f[5], f[2], int(f[3])
            if ctg not in genome: return None
            q = t = sc = alen = clip = 0; t = pos - 1
            for ln, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
                ln = int(ln)
                if op in 'SH': clip += ln; q += ln if op == 'S' else 0
                elif op == 'I': sc -= OPEN + EXT * (ln - 1); q += ln
                elif op == 'D': sc -= OPEN + EXT * (ln - 1); t += ln
                elif op == 'N': t += ln          # an intron is not a gap
                elif op in 'M=X':
                    g = genome[ctg][t:t+ln].upper(); s2 = seq[q:q+ln].upper()
                    for x, y in zip(g, s2): sc += MATCH if x == y else MIS
                    alen += ln; q += ln; t += ln
            return sc, alen, clip
        rows = [(k, score(r[k][0]), score(p[k][0])) for k in cand]
        rows = [(k, a, b) for k, a, b in rows if a and b]
        if rows:
            rs = [a[0] for _, a, b in rows]; ps = [b[0] for _, a, b in rows]
            hi = sum(1 for a, b in zip(rs, ps) if b > a)
            lo = sum(1 for a, b in zip(rs, ps) if b < a)
            print(f"  alignment score on the {len(rows)} reads aligned differently:")
            print(f"     total   reference {sum(rs):10d}   port {sum(ps):10d}"
                  f"   ({100*(sum(ps)-sum(rs))/abs(sum(rs)):+.2f}%)")
            print(f"     port higher on {hi}, lower on {lo}, equal on {len(rows)-hi-lo}")
            print(f"     aligned bases  reference {sum(a[1] for _, a, b in rows):9d}"
                  f"   port {sum(b[1] for _, a, b in rows):9d}")
            print(f"     clipped bases  reference {sum(a[2] for _, a, b in rows):9d}"
                  f"   port {sum(b[2] for _, a, b in rows):9d}")
            print()
