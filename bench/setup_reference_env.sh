#!/usr/bin/env bash
# Create the pinned reference environment the equivalence harness measures against.
#
# The reference is the PYTHON implementation. Everything the harness records as a
# golden is produced by this environment and nothing else, so it is pinned and
# scripted rather than described in prose.
#
# Usage:  bench/setup_reference_env.sh [env_name]
set -euo pipefail

ENV_NAME="${1:-ultra_ref}"

command -v conda >/dev/null 2>&1 || {
  echo "conda not found on PATH. Install miniforge3 first." >&2; exit 1; }

# ---------------------------------------------------------------------------
# Why these packages, and why python 3.12
#
#   parasail-python  4 call sites in modules/help_functions.py
#   python-edlib     2 call sites (classify_read_with_mams.py, help_functions.py)
#   pysam            reads minimap2's SAM, writes indexed.sam / unindexed.sam
#   dill             the 20 on-disk index pickles (see PORTING.md, Hypothesis 3)
#   intervaltree     modules/prefilter_genomic_reads.py
#   gffutils         the GTF -> sqlite index
#   minimap2         live subprocess
#
# python 3.12: 3.12 changed sum() over floats to a compensated summation, so an
# older interpreter is a DIFFERENT reference. Pin it rather than discover it.
# ---------------------------------------------------------------------------
echo "==> creating conda env '$ENV_NAME'"
conda create --yes -n "$ENV_NAME" -c conda-forge -c bioconda \
  python=3.12 pip \
  parasail-python python-edlib pysam dill intervaltree gffutils \
  minimap2

# ---------------------------------------------------------------------------
# namfinder is NOT installable from bioconda on osx-arm64.
#   linux-64 / linux-aarch64 / osx-64 : namfinder 0.1.3 is published
#   osx-arm64                         : absent  -> ksahlin/namfinder#1
# It builds cleanly from source there, so do that rather than fail.
# ---------------------------------------------------------------------------
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

if conda install --yes -c conda-forge -c bioconda namfinder 2>/dev/null; then
  echo "==> namfinder installed from bioconda"
else
  echo "==> namfinder not available for this subdir; building v0.1.3 from source"
  BUILD_DIR="$(mktemp -d)"
  git clone --depth 1 --branch v0.1.3 https://github.com/ksahlin/namfinder "$BUILD_DIR/namfinder"
  # NOTE: no -march=native. INSTALL.sh uses it and that makes the binary
  # non-portable; see PORTING.md Finding 1.
  cmake -B "$BUILD_DIR/namfinder/build" -S "$BUILD_DIR/namfinder" -DCMAKE_BUILD_TYPE=Release
  make -j"$(getconf _NPROCESSORS_ONLN)" -C "$BUILD_DIR/namfinder/build"
  # cmake writes build/namfinder, NOT ./namfinder -- this is Finding 1.
  install -m 0755 "$BUILD_DIR/namfinder/build/namfinder" "$CONDA_PREFIX/bin/namfinder"
  rm -rf "$BUILD_DIR"
fi

echo
echo "==> resolved versions"
python - <<'PY'
import sys
print(f"  python           {sys.version.split()[0]}")
for mod, label in [('parasail','parasail'), ('edlib','edlib'), ('pysam','pysam'),
                   ('dill','dill'), ('intervaltree','intervaltree'), ('gffutils','gffutils')]:
    try:
        m = __import__(mod)
        print(f"  {label:16s} {getattr(m, '__version__', '(no __version__)')}")
    except Exception as e:
        print(f"  {label:16s} MISSING: {e}")
PY
printf "  minimap2         %s\n" "$(minimap2 --version)"
printf "  namfinder        %s\n" "$(namfinder --version 2>&1 | head -1)"

cat <<EOM

==> done.

    conda activate $ENV_NAME
    export REF_PYTHON="\$(conda info --base)/envs/$ENV_NAME/bin/python"

    Every reference run MUST set PYTHONHASHSEED=0; equivalence.sh does this for
    you. Without it 14 of the index outputs differ run to run. See PORTING.md,
    "Determinism".
EOM
