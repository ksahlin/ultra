#!/usr/bin/env bash
# Rewrite the repository's history to drop the committed datasets.
#
#   packaging/history-rewrite.sh [outdir]
#
# Works on a MIRROR CLONE in outdir (default packaging/.rewrite). It never
# touches the repository you run it from and it never pushes. The push is a
# separate, deliberate act -- the commands are printed at the end.
#
# WHAT IT REMOVES: every blob that has only ever lived under data/ and is not
# in the tree of any branch tip. That is 204 blobs, 2 213.6 MB.
#
# WHAT IT KEEPS:
#   - every commit, including the 78 that become empty (--prune-empty never),
#     with author, email, both dates and message unchanged
#   - every branch and tag
#   - the 26 data/ files that exist at the branch tips, which is what
#     README.md links to (data/images)
#   - every non-data/ blob ever committed, including the `torkel` and
#     `torkel2` dev scripts deleted in 2020
#
# THE ONE THING THAT CHANGES CONTENT: tag v0.0.1, whose tree carried 34 data
# files totalling 1 050.6 MB -- nearly half the repository. Its non-data tree
# is unchanged. The other four tags lose nothing at all, because their data/
# files are the same blobs the branch tips still hold.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/packaging/.rewrite}"
MIRROR="$OUT/ultra.git"

command -v git-filter-repo >/dev/null 2>&1 || {
  echo "git-filter-repo is not installed: pip install git-filter-repo" >&2; exit 1; }

rm -rf "$OUT"; mkdir -p "$OUT"
echo "== mirroring $ROOT"
git clone --mirror --no-local "$ROOT" "$MIRROR" 2>&1 | tail -1

# --- record the before state, so the rewrite can be checked against it ------
for r in master develop; do
  git -C "$MIRROR" log --format='%an%x09%ae%x09%aI%x09%cI%x09%s' "$r" > "$OUT/before-$r.txt"
  git -C "$MIRROR" rev-parse "$r^{tree}" > "$OUT/tree-$r.txt"
done
BEFORE=$(git -C "$MIRROR" count-objects -vH | awk '/size-pack/{print $2" "$3}')

# --- build the strip list ---------------------------------------------------
# Keep anything reachable from a BRANCH tip. Deliberately not "any ref tip":
# tag v0.0.1 holds 1 050.6 MB of data/, so including tags in the keep set
# quietly halves the saving, which is exactly what happened on the first run.
echo "== computing the strip list"
git -C "$MIRROR" rev-list --objects --all > "$OUT/objs.txt"
git -C "$MIRROR" cat-file --batch-all-objects \
    --batch-check='%(objectname) %(objecttype) %(objectsize)' > "$OUT/sizes.txt"

MIRROR="$MIRROR" OUT="$OUT" python3 - <<'PY'
import collections, os, subprocess
M, OUT = os.environ['MIRROR'], os.environ['OUT']
size = {}
for line in open(f"{OUT}/sizes.txt"):
    o, t, s = line.split()
    if t == 'blob': size[o] = int(s)
paths = collections.defaultdict(set)
for line in open(f"{OUT}/objs.txt"):
    p = line.rstrip('\n').split(' ', 1)
    if len(p) == 2 and p[1] and p[0] in size: paths[p[0]].add(p[1])

refs = subprocess.run(['git', '-C', M, 'for-each-ref', '--format=%(refname)'],
                      capture_output=True, text=True).stdout.split()
keep = set()
for r in (x for x in refs if not x.startswith('refs/tags/')):
    out = subprocess.run(['git', '-C', M, 'ls-tree', '-r', r],
                         capture_output=True, text=True).stdout
    for l in out.splitlines():
        f = l.split()
        if len(f) >= 4 and f[1] == 'blob': keep.add(f[2])

strip = sorted(o for o, ps in paths.items()
               if o not in keep and all(p.startswith('data/') for p in ps))
# refuse to run on a list that could touch anything outside data/
assert not (set(strip) & keep), "a branch-reachable blob is in the strip list"
assert all(all(p.startswith('data/') for p in paths[o]) for o in strip)
open(f"{OUT}/strip-blobs.txt", "w").write("\n".join(strip) + "\n")
print(f"   {len(strip)} blobs, {sum(size[o] for o in strip)/2**20:.1f} MB")
PY

# --- rewrite ----------------------------------------------------------------
echo "== rewriting"
( cd "$MIRROR" && git filter-repo --strip-blobs-with-ids "$OUT/strip-blobs.txt" \
    --prune-empty never --force >/dev/null )

# --- verify -----------------------------------------------------------------
echo "== verifying"
rc=0
for r in master develop; do
  git -C "$MIRROR" log --format='%an%x09%ae%x09%aI%x09%cI%x09%s' "$r" > "$OUT/after-$r.txt"
  if diff -q "$OUT/before-$r.txt" "$OUT/after-$r.txt" >/dev/null; then
    echo "   $r: all $(wc -l < "$OUT/after-$r.txt" | tr -d ' ') commits, metadata unchanged"
  else echo "   $r: COMMIT METADATA CHANGED" >&2; rc=1; fi
  if [ "$(git -C "$MIRROR" rev-parse "$r^{tree}")" = "$(cat "$OUT/tree-$r.txt")" ]; then
    echo "   $r: HEAD tree identical"
  else echo "   $r: HEAD TREE CHANGED" >&2; rc=1; fi
done
git -C "$MIRROR" cat-file --batch-all-objects --batch-check='%(objectname)' | sort > "$OUT/after-objs.txt"
left=$(comm -12 <(sort "$OUT/strip-blobs.txt") "$OUT/after-objs.txt" | wc -l | tr -d ' ')
[ "$left" = 0 ] && echo "   no stripped blob survives" || { echo "   $left STRIPPED BLOBS SURVIVE" >&2; rc=1; }
[ "$rc" = 0 ] || { echo "VERIFICATION FAILED" >&2; exit 1; }

AFTER=$(git -C "$MIRROR" count-objects -vH | awk '/size-pack/{print $2" "$3}')
cat <<EOF

  pack: $BEFORE  ->  $AFTER
  rewritten mirror: $MIRROR

Nothing has been pushed. To publish it -- which REWRITES PUBLIC HISTORY and
breaks every existing clone -- review it first, then:

    git -C $MIRROR push --force --mirror https://github.com/ksahlin/ultra.git

Before doing that: tell collaborators, and note that GitHub keeps the old
objects reachable by SHA until its garbage collection runs, so anything
genuinely secret in the old history is not destroyed by this alone.
EOF
