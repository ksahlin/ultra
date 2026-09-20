#!/usr/bin/env bash
# Resolve a corpus row from bench/corpora.tsv to local paths, verifying every
# digest first.
#
#   bench/corpora_resolve.sh <name>      -> prints "ref<TAB>gtf<TAB>reads"
#   bench/corpora_resolve.sh --list      -> prints "name<TAB>tier"
#   bench/corpora_resolve.sh --check-all -> verifies every row it can resolve
#
# WHY THIS EXISTS. Rule 3 of this port: check the oracle points at the right
# reference. A harness that silently runs against a substituted or truncated
# corpus is green for the wrong reason, and you do not find out for months.
# So nothing here is optional: if a digest does not match, we exit non-zero and
# say which file and what we got.
#
# $ULTRA_DATA        where the public corpora live      (default $HOME/data)
# $ULTRA_BENCH_CACHE where fetched URLs are cached      (default $HOME/.cache/ultra-bench)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TSV="$ROOT/bench/corpora.tsv"
ULTRA_DATA="${ULTRA_DATA:-$HOME/data}"
CACHE="${ULTRA_BENCH_CACHE:-$HOME/.cache/ultra-bench}"

if command -v shasum >/dev/null 2>&1; then SHA() { shasum -a 256 "$1" | cut -d' ' -f1; }
elif command -v sha256sum >/dev/null 2>&1; then SHA() { sha256sum "$1" | cut -d' ' -f1; }
else echo "need shasum or sha256sum" >&2; exit 1; fi

die() { echo "corpora_resolve: $*" >&2; exit 1; }

# Read the row for a corpus into the global R_* variables.
read_row() {
  local want="$1" found=0
  while IFS=$'\t' read -r name tier ref refsha gtf gtfsha reads readsha notes; do
    case "$name" in ''|'#'*|name) continue ;; esac
    [[ "$name" == "$want" ]] || continue
    R_NAME="$name"; R_TIER="$tier"
    R_REF="$ref";   R_REFSHA="$refsha"
    R_GTF="$gtf";   R_GTFSHA="$gtfsha"
    R_READS="$reads"; R_READSHA="$readsha"
    found=1; break
  done < "$TSV"
  [[ $found -eq 1 ]] || die "no such corpus: $want (see --list)"
}

# Turn a source field into a local path, fetching and caching URLs.
# A .gz URL is decompressed; the digest is checked against the file AS FETCHED.
materialise() {
  local src="$1" sha="$2" kind="$3" out
  case "$src" in
    http://*|https://*)
      mkdir -p "$CACHE"
      local base; base="$(basename "$src")"
      local dl="$CACHE/$base"
      if [[ ! -f "$dl" ]]; then
        echo "corpora_resolve: fetching $kind $src" >&2
        curl -fsSL -o "$dl.part" "$src" || die "fetch failed: $src"
        mv "$dl.part" "$dl"
      fi
      local got; got="$(SHA "$dl")"
      [[ "$got" == "$sha" ]] || die "digest mismatch for $kind $dl
  expected $sha
  got      $got
  (delete the cached file to re-fetch)"
      if [[ "$dl" == *.gz ]]; then
        out="${dl%.gz}"
        [[ -f "$out" ]] || gzip -dc "$dl" > "$out.part" && mv -f "$out.part" "$out" 2>/dev/null || true
        [[ -f "$out" ]] || die "could not decompress $dl"
      else
        out="$dl"
      fi
      ;;
    -)  echo "-"; return 0 ;;
    *)
      out="${src/\$ULTRA_DATA/$ULTRA_DATA}"
      [[ "$out" = /* ]] || out="$ROOT/$out"
      [[ -f "$out" ]] || die "$kind not found: $out
  (set ULTRA_DATA if your public corpora live elsewhere; it is currently $ULTRA_DATA)"
      local got; got="$(SHA "$out")"
      [[ "$got" == "$sha" ]] || die "digest mismatch for $kind $out
  expected $sha
  got      $got"
      ;;
  esac
  echo "$out"
}

case "${1:-}" in
  --list)
    awk -F'\t' '!/^#/ && NF>3 && $1!="name" {printf "%s\t%s\n", $1, $2}' "$TSV"
    exit 0 ;;
  --check-all)
    rc=0
    while IFS=$'\t' read -r name tier _rest; do
      case "$name" in ''|'#'*|name) continue ;; esac
      if out="$("$0" "$name" 2>&1)"; then
        printf 'ok      %-16s %s\n' "$name" "$tier"
      else
        printf 'FAIL    %-16s %s\n' "$name" "$tier"; echo "$out" | sed 's/^/        /'
        rc=1
      fi
    done < "$TSV"
    exit $rc ;;
  ''|-h|--help)
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
esac

read_row "$1"
ref="$(materialise "$R_REF"   "$R_REFSHA"  reference)"
gtf="$(materialise "$R_GTF"   "$R_GTFSHA"  annotation)"
reads="$(materialise "$R_READS" "$R_READSHA" reads)"
printf '%s\t%s\t%s\n' "$ref" "$gtf" "$reads"
