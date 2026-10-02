#!/usr/bin/env bash
# Build the portable uLTRA binaries for every supported target.
#
#   packaging/build-release.sh [outdir]
#
# Produces one .tar.gz per target plus SHA256SUMS in <outdir> (default
# packaging/dist). It does NOT tag, push or publish anything.
#
# WHY .tar.gz AND NOT .zip -- this is not a style choice, see PORTING.md
# Finding 42. A macOS binary that arrives with the com.apple.quarantine
# attribute is killed outright (SIGKILL, exit 137), unsigned and unnotarized
# as these are. Archive Utility, which is what double-clicking a .zip runs,
# PROPAGATES that attribute to the extracted files. Command-line `tar` does
# not. Shipping .zip would ship something that cannot run.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/packaging/dist}"
VERSION="$(grep -m1 '^version' "$ROOT/rust/Cargo.toml" | sed 's/.*"\(.*\)".*/\1/')"

# glibc 2.17 is RHEL/CentOS 7. Older than anything still in service, and the
# floor manylinux2014 uses, so the binary runs on essentially any live Linux.
TARGETS=(
  "x86_64-unknown-linux-gnu.2.17|x86_64-unknown-linux-gnu|linux-x86_64"
  "aarch64-unknown-linux-gnu.2.17|aarch64-unknown-linux-gnu|linux-aarch64"
  "x86_64-apple-darwin|x86_64-apple-darwin|macos-x86_64"
  "aarch64-apple-darwin|aarch64-apple-darwin|macos-arm64"
)

say() { printf '\n== %s\n' "$*"; }

# ---------------------------------------------------------------------------
# zig, for the Linux targets. Installed into a throwaway venv rather than
# expected on the system, so this script needs nothing but rustup and python.
# ---------------------------------------------------------------------------
need_zig() {
  if command -v zig >/dev/null 2>&1; then return; fi
  local venv="$ROOT/packaging/.zigenv"
  [[ -x "$venv/bin/python" ]] || python3 -m venv "$venv"
  "$venv/bin/python" -c 'import ziglang' 2>/dev/null || "$venv/bin/pip" install -q ziglang
  export PATH="$venv/bin:$PATH"
}

command -v cargo-zigbuild >/dev/null 2>&1 || {
  echo "cargo-zigbuild is not installed: cargo install cargo-zigbuild" >&2; exit 1; }
need_zig

mkdir -p "$OUT"
rm -f "$OUT"/*.tar.gz "$OUT"/SHA256SUMS

for spec in "${TARGETS[@]}"; do
  IFS='|' read -r build_target rust_target label <<< "$spec"
  say "$label"
  if [[ "$label" == macos-* ]]; then
    [[ "$(uname -s)" == Darwin ]] || { echo "   skipped: needs a macOS host"; continue; }
    ( cd "$ROOT/rust" && cargo build --release --locked --target "$rust_target" )
  else
    ( cd "$ROOT/rust" && cargo zigbuild --release --locked --target "$build_target" )
  fi

  bin="$ROOT/rust/target/$rust_target/release/uLTRA"
  [[ -f "$bin" ]] || { echo "   no binary at $bin" >&2; exit 1; }

  # Check the glibc floor rather than trusting the target triple to have been
  # honoured. A binary that silently needs a newer glibc fails on the user's
  # machine and nowhere else.
  if [[ "$label" == linux-* ]] && command -v objdump >/dev/null 2>&1; then
    hi="$(objdump -T "$bin" | grep -oE 'GLIBC_[0-9]+\.[0-9]+' | sort -u -t. -k2,2n | tail -1)"
    echo "   highest glibc symbol required: ${hi:-none}"
    [[ "$hi" == "GLIBC_2.17" || -z "$hi" ]] || { echo "   ABOVE THE 2.17 FLOOR" >&2; exit 1; }
  fi

  stage="$(mktemp -d)"
  cp "$bin" "$stage/uLTRA"
  cp "$ROOT/README.md" "$stage/" 2>/dev/null || true
  tar -C "$stage" -czf "$OUT/uLTRA-$VERSION-$label.tar.gz" uLTRA README.md
  rm -rf "$stage"
  echo "   -> $OUT/uLTRA-$VERSION-$label.tar.gz"
done

( cd "$OUT" && shasum -a 256 *.tar.gz > SHA256SUMS )
say "done"
ls -lh "$OUT"
echo
echo "Nothing has been tagged, pushed or published."
