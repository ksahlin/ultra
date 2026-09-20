#!/usr/bin/env bash
# Unit test for equivalence.sh's expected_diff().
#
# WHY THIS IS ITS OWN TEST. expected_diff decides whether a differing structure
# is an approved divergence or a bug. It is the one piece of the harness whose
# failure mode is SILENT: if it wrongly returns "expected", a real regression is
# printed as DIVERGED and scrolls past as though it were fine. And it cannot be
# exercised through `stage index verify` until the port implements dump-index,
# so it is tested directly.
#
# The load-bearing case is the second one: the Finding 22 divergence is approved
# on droso-20k ONLY. On SIRV the port must still match, because every SIRV part
# has exactly one active gene - that is what makes SIRV the control.
set -u
BENCH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

expected_diff() {
  awk -F'\t' -v st="$1" -v sr="$2" -v co="$3" '
    !/^#/ && NF>3 && $1!="stage" && $1==st && $2==sr {
      if ($4=="*") { found=1; exit }
      n=split($4, a, ","); for (i=1;i<=n;i++) if (a[i]==co) { found=1; exit }
    }
    END { exit(found?0:1) }' "$BENCH/stage_diffs.tsv"
}

fails=0
t() {
  local desc="$1" expect="$2"; shift 2
  local got; if expected_diff "$@"; then got=yes; else got=no; fi
  if [[ "$got" == "$expect" ]]; then
    printf '  ok    %s\n' "$desc"
  else
    printf '  FAIL  %s  (expected=%s got=%s)\n' "$desc" "$expect" "$got"; fails=$((fails+1))
  fi
}

t "the listed divergence, on its own corpus"       yes index gene_to_small_segments droso-20k
t "SAME structure on sirv-10k is NOT excused"      no  index gene_to_small_segments sirv-10k
t "SAME structure on smoke is NOT excused"         no  index gene_to_small_segments smoke
t "a different structure on droso is NOT excused"  no  index segment_to_gene        droso-20k
t "an unknown structure is NOT excused"            no  index nonesuch               droso-20k
t "an unknown stage is NOT excused"                no  align gene_to_small_segments droso-20k

echo
if [[ $fails -eq 0 ]]; then echo "stage_diffs: all checks passed"; else echo "stage_diffs: $fails FAILED"; fi
exit $fails
