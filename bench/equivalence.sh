#!/usr/bin/env bash
# uLTRA equivalence harness.
#
#   bench/equivalence.sh check              environment + oracle sanity (RUN THIS FIRST)
#   bench/equivalence.sh record [corpus]    run the REFERENCE, write golden manifests
#   bench/equivalence.sh verify [corpus]    run the PORT, compare against the goldens
#   bench/equivalence.sh stable [corpus]    record twice and diff - proves the goldens
#                                           are reproducible before you trust them
#   bench/equivalence.sh cli record|verify  the CLI contract: exit codes, stdout,
#                                           stderr and the outfolder side effect
#   bench/equivalence.sh stage index record|verify [corpus]
#                                           the index stage oracle - a semantic
#                                           contract, because the reference's
#                                           index is Python pickles + sqlite
#   bench/equivalence.sh list               show cases and corpora
#
# Environment:
#   REF_PYTHON          interpreter with parasail, edlib, pysam, dill, gffutils
#                       (default: $HOME/miniforge3/envs/ultra_ref/bin/python)
#   PORT_BIN            the Rust binary          (default: rust/target/release/ultra)
#   ULTRA_DATA          public corpora root      (default: $HOME/data)
#   CORPORA             tier to run: smoke|core|heavy   (default: core)
#   CASES               regex selecting case names      (default: all)
#
# ---------------------------------------------------------------------------
# THE FOUR RULES THIS SCRIPT ENCODES. Each cost days somewhere else.
#
# 1. Never verify through a privileged path. Not this script's job, but the
#    release checks must download through a browser, not `gh release download`.
# 2. Never derive a number by reimplementing the tool's logic. Everything here
#    is read off an actual run.
# 3. Check the oracle points at the right reference. `check` refuses to run if
#    REF_PYTHON is missing, cannot import the dependencies, or is not the
#    interpreter you think it is - and corpora_resolve.sh verifies every input
#    digest, so a substituted corpus cannot make the suite green.
# 4. Guard against stale binaries. check_bin_fresh HARD-FAILS if any Rust source
#    is newer than PORT_BIN.
# ---------------------------------------------------------------------------
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH="$ROOT/bench"
GOLDEN="$BENCH/golden"
REF_PYTHON="${REF_PYTHON:-$HOME/miniforge3/envs/ultra_ref/bin/python}"
PORT_BIN="${PORT_BIN:-$ROOT/rust/target/release/uLTRA}"
CORPORA="${CORPORA:-core}"
CASES="${CASES:-.}"
WORK="${ULTRA_BENCH_WORK:-${TMPDIR:-/tmp}/ultra-equiv.$$}"

if command -v shasum >/dev/null 2>&1; then SHA() { shasum -a 256 "$1" | cut -d' ' -f1; }
elif command -v sha256sum >/dev/null 2>&1; then SHA() { sha256sum "$1" | cut -d' ' -f1; }
else echo "need shasum or sha256sum" >&2; exit 1; fi

RED=''; GRN=''; YEL=''; OFF=''
if [[ -t 1 ]]; then RED=$'\033[31m'; GRN=$'\033[32m'; YEL=$'\033[33m'; OFF=$'\033[0m'; fi
ok()   { printf '%sok%s      %s\n'   "$GRN" "$OFF" "$*"; }
bad()  { printf '%sFAIL%s    %s\n'   "$RED" "$OFF" "$*"; }
warn() { printf '%swarn%s    %s\n'   "$YEL" "$OFF" "$*"; }
die()  { echo "equivalence: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Output classification. This is the contract, and it is the most important
# thing in the file. See PORTING.md Part 5.
#
#   hash      byte-for-byte.
#   sam       header byte-for-byte, alignment records SORTED before hashing.
#             uLTRA emits records in worker-completion order, which is
#             scheduling-dependent: measured, four concurrent runs of the same
#             input at the default --t 3 produced THREE different files, all
#             with identical content in a different order. At --t 1 it is
#             stable. So record order is not a contract the reference offers,
#             and hashing it would make the suite fail on a busy machine.
#             See PORTING.md Finding 32.
#   gunzip    hash the DECOMPRESSED bytes. gzip stores an mtime in header bytes
#             5-8, so seeds.txt.gz and reads_tmp.fa.gz differ every run while
#             their contents are identical. Measured, not assumed.
#   scrub_pg  hash with @PG and @HD removed. minimap2's SAMs carry its version
#             and the verbatim command line, including absolute paths.
#   exists    record presence and size only. Logs whose content is timings.
#   oracle    NOT hashed. Python pickles and the gffutils sqlite database: a
#             Rust port will never produce these, so byte-identity is not even
#             a coherent goal. Covered by dump_reference.py stage oracles.
# ---------------------------------------------------------------------------
classify() {
  case "$1" in
    *.pickle|database.db)                      echo oracle ;;
    *.gz)                                      echo gunzip ;;
    minimap2.sam|indexed.sam|unindexed.sam)    echo scrub_pg ;;
    *.sam)                                     echo sam ;;
    minimap2_errors.*|namfinder_stderr.*|*.stderr) echo exists ;;
    *)                                         echo hash ;;
  esac
}

digest_for() {
  local f="$1" kind="$2"
  case "$kind" in
    hash)     SHA "$f" ;;
    gunzip)   gzip -dc "$f" 2>/dev/null | { if command -v shasum >/dev/null 2>&1; then shasum -a 256; else sha256sum; fi; } | cut -d' ' -f1 ;;
    scrub_pg) grep -v '^@PG' "$f" | grep -v '^@HD' | { if command -v shasum >/dev/null 2>&1; then shasum -a 256; else sha256sum; fi; } | cut -d' ' -f1 ;;
    sam)      { grep '^@' "$f"; grep -v '^@' "$f" | LC_ALL=C sort; } | { if command -v shasum >/dev/null 2>&1; then shasum -a 256; else sha256sum; fi; } | cut -d' ' -f1 ;;
    exists)   echo "size:$(wc -c < "$f" | tr -d ' ')" ;;
    oracle)   echo "oracle" ;;
  esac
}

manifest_of() {   # manifest_of <dir>  -> "kind<TAB>digest<TAB>relpath" sorted
  local d="$1"
  ( cd "$d" && find . -type f | sed 's|^\./||' | LC_ALL=C sort ) | while read -r rel; do
    local kind; kind="$(classify "$(basename "$rel")")"
    printf '%s\t%s\t%s\n' "$kind" "$(digest_for "$d/$rel" "$kind")" "$rel"
  done
}

# ---------------------------------------------------------------------------
corpora_of_tier() {
  local want="$1"
  "$BENCH/corpora_resolve.sh" --list | while IFS=$'\t' read -r name tier; do
    case "$want:$tier" in
      smoke:smoke|core:smoke|core:core|heavy:*) echo "$name" ;;
    esac
  done
}

each_case() {   # emits: name<TAB>entry<TAB>tier<TAB>needs<TAB>args
  awk -F'\t' '!/^#/ && NF>3 && $1!="name" {printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5}' "$BENCH/cases.tsv"
}

case_applies() {   # case_tier corpus_tier -> 0 if the case should run
  local ct="$1" pt="$2"
  case "$ct:$pt" in
    smoke:*)        return 0 ;;
    core:core|core:heavy) return 0 ;;
    heavy:heavy)    return 0 ;;
    *)              return 1 ;;
  esac
}

tier_of() { "$BENCH/corpora_resolve.sh" --list | awk -F'\t' -v n="$1" '$1==n{print $2}'; }

# ---------------------------------------------------------------------------
check_bin_fresh() {   # RULE 4
  [[ -x "$PORT_BIN" ]] || die "PORT_BIN not executable: $PORT_BIN
  build it first, or set PORT_BIN."
  # RULE 3, the case-sensitivity trap. macOS resolves target/release/ultra to
  # the real uLTRA, so a wrong default passes here and fails on Linux. Compare
  # the basename against the directory listing, which IS case-sensitive.
  local d b; d="$(dirname "$PORT_BIN")"; b="$(basename "$PORT_BIN")"
  if ! ls -1 "$d" 2>/dev/null | grep -qx -- "$b"; then
    die "PORT_BIN resolves but its name does not match on a case-sensitive
  filesystem: $PORT_BIN
  the directory contains: $(ls -1 "$d" | grep -ix -- "$b" | tr '\n' ' ')"
  fi
  # Only meaningful when PORT_BIN really is built from rust/src. The harness's
  # own test doubles (bench's deliberately-broken "ports") live elsewhere and
  # are not compiled from anything, so comparing mtimes against them is a false
  # alarm -- one that fired during development and cost a confusing few minutes.
  case "$PORT_BIN" in
    "$ROOT/rust/target/"*) ;;
    *) return 0 ;;
  esac
  local newer
  newer="$(find "$ROOT/rust/src" "$ROOT/rust/Cargo.toml" -newer "$PORT_BIN" 2>/dev/null | head -5)"
  if [[ -n "$newer" ]]; then
    bad "STALE BINARY: these are newer than $PORT_BIN"
    echo "$newer" | sed 's/^/        /'
    die "refusing to verify against a stale build. cargo build --release first."
  fi
}

# ---------------------------------------------------------------------------
cmd_check() {   # RULE 3
  local rc=0
  echo "== reference interpreter =="
  if [[ -x "$REF_PYTHON" ]]; then ok "REF_PYTHON $REF_PYTHON"; else bad "REF_PYTHON not executable: $REF_PYTHON"; rc=1; fi
  if [[ -x "$REF_PYTHON" ]]; then
    "$REF_PYTHON" - <<'PY' || rc=1
import sys
v = sys.version_info
print(f"        python {v.major}.{v.minor}.{v.micro}")
if (v.major, v.minor) < (3, 12):
    print("        FAIL: python < 3.12 is a DIFFERENT reference (sum() over floats changed)")
    raise SystemExit(1)
missing = []
for m in ("parasail", "edlib", "pysam", "dill", "intervaltree", "gffutils"):
    try: __import__(m)
    except Exception as e: missing.append(f"{m} ({e})")
if missing:
    print("        FAIL missing:", ", ".join(missing)); raise SystemExit(1)
print("        parasail edlib pysam dill intervaltree gffutils all import")
PY
  fi
  echo "== the reference is THIS repository, not an installed copy =="
  if [[ -f "$ROOT/uLTRA" && -d "$ROOT/modules" ]]; then ok "$ROOT/uLTRA + modules/"; else bad "no uLTRA script at $ROOT"; rc=1; fi
  if command -v uLTRA >/dev/null 2>&1; then
    warn "an installed uLTRA is also on PATH ($(command -v uLTRA)); this harness ignores it"
  fi
  echo "== external tools =="
  for t in minimap2 namfinder; do
    if command -v "$t" >/dev/null 2>&1; then ok "$t $("$t" --version 2>&1 | head -1)"; else bad "$t not on PATH"; rc=1; fi
  done
  echo "== corpora =="
  "$BENCH/corpora_resolve.sh" --check-all || rc=1
  echo "== cases =="
  ok "$(each_case | wc -l | tr -d ' ') cases in bench/cases.tsv"
  return $rc
}

# ---------------------------------------------------------------------------
# run_case <engine> <corpus> <case_name> <entry> <needs> <args> <outdir>
# engine is either "ref" or "port".
run_case() {
  local engine="$1" corpus="$2" cname="$3" entry="$4" needs="$5" args="$6" out="$7"
  local resolved ref gtf reads
  resolved="$("$BENCH/corpora_resolve.sh" "$corpus")" || return 90
  IFS=$'\t' read -r ref gtf reads <<< "$resolved"

  mkdir -p "$out"
  local -a cmd
  if [[ "$engine" == ref ]]; then cmd=(env PYTHONHASHSEED=0 "$REF_PYTHON" "$ROOT/uLTRA")
  else                            cmd=("$PORT_BIN"); fi

  local idxdir="$WORK/idx/$corpus/$needs"
  case "$entry" in
    index)    ( cd "$ROOT" && "${cmd[@]}" index $args "$ref" "$gtf" "$out" ) ;;
    align)    ( cd "$ROOT" && "${cmd[@]}" align $args --index "$idxdir" "$ref" "$reads" "$out" ) ;;
    pipeline) ( cd "$ROOT" && "${cmd[@]}" pipeline $args "$ref" "$gtf" "$reads" "$out" ) ;;
    *) die "unknown entry: $entry" ;;
  esac >"$out.stdout" 2>"$out.stderr"
  echo $? > "$out.exit"
}

# Build the index a set of align cases depends on, once per corpus.
ensure_index() {
  local engine="$1" corpus="$2" needs="$3"
  [[ "$needs" == "-" ]] && return 0
  local idxdir="$WORK/idx/$corpus/$needs"
  [[ -d "$idxdir" && -f "$idxdir/database.db" ]] && return 0
  local iargs
  iargs="$(each_case | awk -F'\t' -v n="$needs" '$1==n{print $5}')"
  mkdir -p "$idxdir"
  run_case "$engine" "$corpus" "$needs" index - "$iargs" "$idxdir"
  local rc; rc="$(cat "$idxdir.exit" 2>/dev/null || echo 99)"
  [[ "$rc" == 0 ]] || { bad "index case '$needs' failed on $corpus (exit $rc); align cases needing it are skipped"; return 1; }
}

# ---------------------------------------------------------------------------
sweep() {   # sweep <record|verify> [corpus]
  local mode="$1" only="${2:-}"
  local engine=ref; [[ "$mode" == verify ]] && engine=port
  [[ "$engine" == port ]] && check_bin_fresh

  mkdir -p "$WORK"
  local total=0 pass=0 fail=0 skip=0

  local corpora; corpora="$(corpora_of_tier "$CORPORA")"
  [[ -n "$only" ]] && corpora="$only"

  for corpus in $corpora; do
    local ptier; ptier="$(tier_of "$corpus")"
    echo
    echo "######## corpus: $corpus (tier $ptier) ########"
    while IFS=$'\t' read -r cname entry ctier needs args; do
      [[ "$cname" =~ $CASES ]] || continue
      case_applies "$ctier" "$ptier" || continue
      total=$((total+1))

      if [[ "$entry" == align ]]; then
        ensure_index "$engine" "$corpus" "$needs" || { skip=$((skip+1)); continue; }
      fi

      local out="$WORK/$mode/$corpus/$cname"
      rm -rf "$out"; mkdir -p "$out"
      run_case "$engine" "$corpus" "$cname" "$entry" "$needs" "$args" "$out"
      local rc; rc="$(cat "$out.exit")"

      local gdir="$GOLDEN/$corpus/$cname"
      if [[ "$mode" == record ]]; then
        mkdir -p "$gdir"
        manifest_of "$out" > "$gdir/manifest.tsv"
        echo "$rc" > "$gdir/exit"
        ok "$cname  (exit $rc, $(wc -l < "$gdir/manifest.tsv" | tr -d ' ') files)"
        pass=$((pass+1))
      else
        if [[ ! -f "$gdir/manifest.tsv" ]]; then
          warn "$cname  no golden recorded; skipping"; skip=$((skip+1)); continue
        fi
        local grc; grc="$(cat "$gdir/exit")"
        local d; d="$(diff <(cat "$gdir/manifest.tsv") <(manifest_of "$out") || true)"
        if [[ "$rc" == "$grc" && -z "$d" ]]; then
          ok "$cname"; pass=$((pass+1))
        else
          bad "$cname"
          [[ "$rc" == "$grc" ]] || echo "        exit $grc -> $rc"
          [[ -z "$d" ]] || echo "$d" | head -20 | sed 's/^/        /'
          fail=$((fail+1))
        fi
      fi
    done < <(each_case)
  done

  echo
  echo "======== $mode: $pass ok, $fail failed, $skip skipped, of $total ========"
  [[ $fail -eq 0 ]]
}

cmd_stable() {   # record twice into two trees and diff: are the goldens reproducible?
  local a="$WORK/stableA" b="$WORK/stableB" rc=0
  for pass in A B; do
    local dir="$WORK/stable$pass"
    GOLDEN="$dir" sweep record "${1:-}" >/dev/null || true
  done
  echo "== stable: comparing two independent recordings =="
  local n=0 d=0
  while read -r f; do
    n=$((n+1))
    local rel="${f#$a/}"
    if ! diff -q "$f" "$b/$rel" >/dev/null 2>&1; then
      bad "unstable: $rel"; diff "$f" "$b/$rel" | head -8 | sed 's/^/        /'; d=$((d+1))
    fi
  done < <(find "$a" -name manifest.tsv 2>/dev/null)
  echo "======== stable: $((n-d)) of $n manifests reproducible ========"
  [[ $d -eq 0 ]]
}

cmd_list() {
  echo "== corpora =="; "$BENCH/corpora_resolve.sh" --list | sed 's/^/  /'
  echo "== cases =="; each_case | awk -F'\t' '{printf "  %-22s %-9s %-6s %s\n",$1,$2,$3,$5}'
}

# ---------------------------------------------------------------------------
# CLI contract.
#
# Captures exit code, stdout, stderr and the outfolder side effect. Everything
# that legitimately varies between machines and runs is scrubbed; everything
# else is contract.
#
# What is scrubbed and why:
#   the temp outfolder path   -> {OUT}     differs every run
#   the repository root       -> {ROOT}    differs per checkout
#   "line 123," in tracebacks -> "line N," moves whenever the .py file is edited
#   float seconds             -> {T}       timings
# Anything else - wording, argparse usage blocks, the order of the flags in the
# usage line - IS the contract. See PORTING.md Findings 15-18.
# ---------------------------------------------------------------------------
scrub() {   # scrub <outdir>
  local out="$1"
  sed -e "s|$out|{OUT}|g" \
      -e "s|$ROOT|{ROOT}|g" \
      -e "s|$HOME|{HOME}|g" \
      -e 's|line [0-9][0-9]*,|line N,|g' \
      -e 's|[0-9][0-9]*\.[0-9][0-9]*e-[0-9]*|{T}|g' \
      -e 's|[0-9][0-9]*\.[0-9][0-9][0-9][0-9]*|{T}|g'
}

cli_expand() {   # cli_expand <args> <outdir>
  local a="$1" out="$2"
  a="${a//\{OUT\}/$out}"
  a="${a//\{REF\}/$ROOT/test/SIRV_genes.fasta}"
  a="${a//\{GTF\}/$ROOT/test/SIRV_genes_C_170612a.gtf}"
  a="${a//\{READS\}/$ROOT/test/reads.fa}"
  echo "$a"
}

each_cli_case() {   # name<TAB>status<TAB>args
  awk -F'\t' '!/^#/ && NF>2 && $1!="name" {printf "%s\t%s\t%s\n",$1,$2,$3}' "$BENCH/cli_cases.tsv"
}

cli_run_one() {   # cli_run_one <engine> <args> <outdir> <dest>
  local engine="$1" args="$2" out="$3" dest="$4"
  local -a cmd
  if [[ "$engine" == ref ]]; then cmd=(env PYTHONHASHSEED=0 COLUMNS=80 "$REF_PYTHON" "$ROOT/uLTRA")
  else                            cmd=(env COLUMNS=80 "$PORT_BIN"); fi
  mkdir -p "$dest"
  ( cd "$ROOT" && "${cmd[@]}" $args ) >"$dest/stdout.raw" 2>"$dest/stderr.raw"
  echo $? > "$dest/exit"
  scrub "$out" < "$dest/stdout.raw" > "$dest/stdout"
  scrub "$out" < "$dest/stderr.raw" > "$dest/stderr"
  rm -f "$dest/stdout.raw" "$dest/stderr.raw"
  # Finding 17: the outfolder side effect is observable, so it is recorded.
  if [[ -d "$out" ]]; then
    printf 'created=yes entries=%s\n' "$(ls -A "$out" 2>/dev/null | wc -l | tr -d ' ')" > "$dest/outfolder"
  else
    printf 'created=no\n' > "$dest/outfolder"
  fi
}

cmd_cli() {   # cmd_cli <record|verify>
  local mode="${1:-record}"
  local engine=ref; [[ "$mode" == verify ]] && engine=port
  [[ "$engine" == port ]] && check_bin_fresh
  local gbase="$GOLDEN/cli"
  local pass=0 fail=0 total=0

  local pending=0
  while IFS=$'\t' read -r cname cstatus args; do
    [[ "$cname" =~ $CASES ]] || continue
    # A pending case is one the port cannot satisfy until a later stage. It is
    # still RUN and still recorded; on verify it is reported, never hidden.
    if [[ "$mode" == verify && "$cstatus" == pending:* ]]; then
      warn "$cname  (${cstatus}) - not expected to match yet"
      pending=$((pending+1)); continue
    fi
    total=$((total+1))
    local out; out="$(mktemp -d "${TMPDIR:-/tmp}/ultra-cli.XXXXXX")"
    rmdir "$out"                      # the tool must create it, not us
    local exp; exp="$(cli_expand "$args" "$out")"
    local dest="$WORK/cli/$mode/$cname"
    rm -rf "$dest"
    cli_run_one "$engine" "$exp" "$out" "$dest"
    rm -rf "$out"

    if [[ "$mode" == record ]]; then
      mkdir -p "$gbase/$cname"
      cp "$dest/exit" "$dest/stdout" "$dest/stderr" "$dest/outfolder" "$gbase/$cname/"
      ok "$cname  (exit $(cat "$dest/exit"), $(cat "$dest/outfolder"))"
      pass=$((pass+1))
    else
      if [[ ! -d "$gbase/$cname" ]]; then warn "$cname  no golden"; continue; fi
      local d=""
      for f in exit stdout stderr outfolder; do
        if ! diff -q "$gbase/$cname/$f" "$dest/$f" >/dev/null 2>&1; then
          d+=$'\n'"    --- $f ---"$'\n'"$(diff "$gbase/$cname/$f" "$dest/$f" | head -12)"
        fi
      done
      if [[ -z "$d" ]]; then ok "$cname"; pass=$((pass+1))
      else bad "$cname"; echo "$d" | sed 's/^/    /'; fail=$((fail+1)); fi
    fi
  done < <(each_cli_case)

  echo
  if [[ "$mode" == verify ]]; then
    echo "======== cli $mode: $pass ok, $fail failed, $pending pending, of $((total+pending)) ========"
  else
    echo "======== cli $mode: $pass ok, $fail failed, of $total ========"
  fi
  [[ $fail -eq 0 ]]
}

# Every case must be classified. Refuse to run a matrix that has drifted.
cli_audit() {
  local bad=0 n=0
  while IFS=$'\t' read -r cname cstatus _args; do
    n=$((n+1))
    case "$cstatus" in
      contract|pending:*) ;;
      *) bad "unclassified CLI case '$cname' (status='$cstatus')"; bad=$((bad+1)) ;;
    esac
  done < <(each_cli_case)
  local c p
  c=$(each_cli_case | awk -F'\t' '$2=="contract"' | wc -l | tr -d ' ')
  p=$(each_cli_case | awk -F'\t' '$2 ~ /^pending:/' | wc -l | tr -d ' ')
  echo "cli cases: $n total = $c contract + $p pending"
  each_cli_case | awk -F'\t' '$2 ~ /^pending:/ {printf "  pending %-10s %s\n", $2, $1}'
  [[ $bad -eq 0 ]]
}

# ---------------------------------------------------------------------------
# Stage oracles.
#
# The index stage has NO byte-identity contract: its output is 20 Python
# pickles and a gffutils sqlite database, neither of which a Rust port will
# ever produce. So the contract is the canonical text rendering produced by
# bench/dump_reference.py (from the reference) and by `uLTRA dump-index`
# (from the port), compared structure by structure.
#
# Renderings are large -- 265 MB on Drosophila -- so the golden stores a
# sha256 per structure, not the text.
# ---------------------------------------------------------------------------
expected_diff() {   # expected_diff <stage> <structure> <corpus> -> 0 if listed
  awk -F'\t' -v st="$1" -v sr="$2" -v co="$3" '
    !/^#/ && NF>3 && $1!="stage" && $1==st && $2==sr {
      if ($4=="*") { found=1; exit }
      n=split($4, a, ","); for (i=1;i<=n;i++) if (a[i]==co) { found=1; exit }
    }
    END { exit(found?0:1) }' "$BENCH/stage_diffs.tsv"
}

cmd_stage() {   # cmd_stage index record|verify [corpus]
  local stage="${1:-index}" mode="${2:-record}" only="${3:-}"
  [[ "$stage" == index ]] || die "unknown stage: $stage (only 'index' exists so far)"
  local engine=ref; [[ "$mode" == verify ]] && engine=port
  [[ "$engine" == port ]] && check_bin_fresh
  mkdir -p "$WORK"

  local corpora; corpora="$(corpora_of_tier "$CORPORA")"
  [[ -n "$only" ]] && corpora="$only"
  local match=0 diverged=0 failed=0 stale=0

  for corpus in $corpora; do
    echo
    echo "######## stage=$stage corpus=$corpus ########"
    local idx="$WORK/stage/$corpus/idx" ren="$WORK/stage/$corpus/$mode"
    rm -rf "$idx" "$ren"; mkdir -p "$idx" "$ren"

    # The reference must build an index directory to be rendered from; the port
    # renders straight from the inputs.
    if [[ "$engine" == ref ]]; then
      run_case "$engine" "$corpus" idx-default index - "" "$idx"
      local rc; rc="$(cat "$idx.exit")"
      if [[ "$rc" != 0 ]]; then bad "$corpus: index build failed (exit $rc)"; failed=$((failed+1)); continue; fi
    fi

    if [[ "$engine" == ref ]]; then
      "$REF_PYTHON" "$BENCH/dump_reference.py" --index "$idx" --out "$ren" >/dev/null || {
        bad "$corpus: dump_reference.py failed"; failed=$((failed+1)); continue; }
    else
      # The port takes the INPUTS, not an index directory. Stage 2 verifies index
      # CONSTRUCTION; the on-disk index format is a separate concern that arrives
      # when `align` needs to load one. Passing the inputs also means the port is
      # never handed a reference-built index it could accidentally read.
      local resolved ref gtf _reads
      resolved="$("$BENCH/corpora_resolve.sh" "$corpus")" || { failed=$((failed+1)); continue; }
      IFS=$'\t' read -r ref gtf _reads <<< "$resolved"
      "$PORT_BIN" dump-index "$ref" "$gtf" "$ren" >/dev/null || {
        bad "$corpus: port dump-index failed"; failed=$((failed+1)); continue; }
    fi

    local gdir="$GOLDEN/$corpus/stage-$stage"
    if [[ "$mode" == record ]]; then
      mkdir -p "$gdir"
      ( cd "$ren" && for f in *.txt; do printf '%s\t%s\n' "${f%.txt}" "$(SHA "$f")"; done ) \
        | LC_ALL=C sort > "$gdir/structures.tsv"
      ok "$corpus  ($(wc -l < "$gdir/structures.tsv" | tr -d ' ') structures)"
      match=$((match+1))
    else
      [[ -f "$gdir/structures.tsv" ]] || { warn "$corpus: no stage golden"; continue; }
      while IFS=$'\t' read -r name want; do
        local got; got="$(SHA "$ren/$name.txt" 2>/dev/null || echo MISSING)"
        if [[ "$got" == "$want" ]]; then
          if expected_diff "$stage" "$name" "$corpus"; then
            printf '%sSTALE%s   %-34s listed as an expected divergence but MATCHES - remove the row\n' "$YEL" "$OFF" "$name"
            stale=$((stale+1))
          else
            match=$((match+1))
          fi
        elif expected_diff "$stage" "$name" "$corpus"; then
          printf '%sDIVERGED%s %-34s (expected; see bench/stage_diffs.tsv)\n' "$YEL" "$OFF" "$name"
          diverged=$((diverged+1))
        else
          bad "$name  differs and is NOT an expected divergence"
          failed=$((failed+1))
        fi
      done < "$gdir/structures.tsv"
    fi
  done

  echo
  echo "======== stage $stage $mode: $match match, $diverged diverged, $failed FAILED, $stale stale ========"
  [[ $failed -eq 0 && $stale -eq 0 ]]
}

case "${1:-check}" in
  check)  cmd_check ;;
  stage)  shift; cmd_stage "${1:-index}" "${2:-record}" "${3:-}" ;;
  cli)    shift; cmd_cli "${1:-record}" ;;
  cli_audit) cli_audit ;;
  record) shift; sweep record "${1:-}" ;;
  verify) shift; sweep verify "${1:-}" ;;
  stable) shift; cmd_stable "${1:-}" ;;
  list)   cmd_list ;;
  *) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
