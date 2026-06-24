#!/usr/bin/env bash
# ============================================================================
# install_hpc.sh — Set up PyMieDAP in a Python virtual environment on a
#                  Linux HPC cluster (x86_64 or aarch64, SLURM / PBS / LSF).
#
# Usage (run from the repo root):
#   bash install_hpc.sh              # auto-detect Python and gfortran
#   bash install_hpc.sh --python python3.11  # force a specific interpreter
#   bash install_hpc.sh --venv /path/to/venv # custom venv location
#
# What this script does:
#   1. Optionally loads HPC modules (Python, GCC/gfortran)
#   2. Creates a Python virtual environment
#   3. Installs all Python dependencies (numpy < 2.0 is enforced)
#   4. Compiles the five Fortran extension modules with f2py
#   5. Verifies the install with a quick smoke test
# ============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
PYTHON_BIN="python3"          # override with --python
VENV_DIR="$(pwd)/venv"        # override with --venv
VENV_ACTIVATE=""
SKIP_MODULES=0                # set 1 if your cluster has no `module` command

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --python)  PYTHON_BIN="$2"; shift 2 ;;
    --venv)    VENV_DIR="$2";   shift 2 ;;
    --no-modules) SKIP_MODULES=1; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

echo "========================================================"
echo " PyMieDAP HPC installer"
echo " Repo  : $REPO_ROOT"
echo " Python: $PYTHON_BIN"
echo " Venv  : $VENV_DIR"
echo "========================================================"

# ── Step 1: Load HPC modules (edit to match your cluster) ────────────────────
if [[ $SKIP_MODULES -eq 0 ]] && command -v module &>/dev/null; then
  echo ""
  echo "[1/5] Loading HPC modules …"
  # ── EDIT THESE to match your cluster's module names ──────────────────────
  # Common variants are shown; uncomment whichever applies.
  #
  # module load python/3.11          # or python3, Python/3.11.x, etc.
  # module load gcc/12               # provides gfortran
  # module load anaconda3            # alternative: use conda Python
  #
  # If your cluster requires specific module combinations, replace the block
  # above.  Run `module avail python` and `module avail gcc` to discover the
  # exact names on your system.
  #
  # If modules are already loaded (e.g. in your .bashrc) or not needed,
  # run this script with --no-modules.
  echo "  (module loading skipped — edit Section 1 of install_hpc.sh to"
  echo "   load Python and GCC modules for your specific cluster)"
fi

# ── Step 2: Verify prerequisites ──────────────────────────────────────────────
echo ""
echo "[2/5] Checking prerequisites …"

if ! command -v "$PYTHON_BIN" &>/dev/null; then
  echo "  ERROR: '$PYTHON_BIN' not found."
  echo "  Load the Python module first (see Step 1 above) or pass --python."
  exit 1
fi

PY_VERSION=$("$PYTHON_BIN" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_MAJOR=$("$PYTHON_BIN" -c "import sys; print(sys.version_info.major)")
PY_MINOR=$("$PYTHON_BIN" -c "import sys; print(sys.version_info.minor)")
echo "  Python  : $("$PYTHON_BIN" --version)"

if [[ $PY_MAJOR -lt 3 || ($PY_MAJOR -eq 3 && $PY_MINOR -lt 9) ]]; then
  echo "  ERROR: Python >= 3.9 required (found $PY_VERSION)."
  exit 1
fi

if ! command -v gfortran &>/dev/null; then
  echo "  ERROR: gfortran not found. Load a GCC module or install gfortran."
  echo "  On RHEL/CentOS: sudo yum install gcc-gfortran"
  echo "  On Debian/Ubuntu: sudo apt-get install gfortran"
  exit 1
fi
echo "  gfortran: $(gfortran --version | head -1)"

# ── Step 3: Create virtual environment ────────────────────────────────────────
echo ""
echo "[3/5] Creating virtual environment at $VENV_DIR …"

if [[ -d "$VENV_DIR" ]]; then
  echo "  Venv already exists — activating existing one."
else
  "$PYTHON_BIN" -m venv "$VENV_DIR"
  echo "  Created."
fi

source "$VENV_DIR/bin/activate"
VENV_ACTIVATE="$VENV_DIR/bin/activate"
echo "  Activated: $(python --version) at $(which python)"

# Upgrade pip silently
pip install --quiet --upgrade pip setuptools wheel

# ── Step 4: Install Python dependencies ───────────────────────────────────────
echo ""
echo "[4/5] Installing Python dependencies …"
echo "  (numpy is pinned to < 2.0 — the Fortran extensions require numpy 1.x)"

pip install --quiet -r "$REPO_ROOT/requirements.txt"

NUMPY_VER=$(python -c "import numpy; print(numpy.__version__)")
echo "  numpy installed: $NUMPY_VER"

NUMPY_MAJOR=$(python -c "import numpy; print(numpy.__version__.split('.')[0])")
if [[ "$NUMPY_MAJOR" == "2" ]]; then
  echo "  ERROR: numpy 2.x was installed despite the pin. The Fortran modules"
  echo "  will not load. Check for conflicting package constraints."
  exit 1
fi

# ── Step 5: Compile Fortran extension modules ─────────────────────────────────
echo ""
echo "[5/5] Compiling Fortran extensions with f2py …"
echo "  Using gfortran: $(which gfortran)"
echo "  Build target : $(python -c 'import sysconfig; print(sysconfig.get_platform())')"
cd "$REPO_ROOT"

compile_module() {
  local name="$1"
  shift
  echo "  Compiling $name …"
  python -m numpy.f2py -c -m "$name" "$@" --f77flags="-O2" --quiet 2>&1 \
    | grep -v "^$\|UserWarning\|FutureWarning\|deprecated" || true
  if python -c "import $name" 2>/dev/null; then
    echo "    OK: $name imported successfully."
  else
    echo "    FAILED: $name could not be imported after compilation."
    exit 1
  fi
}

compile_module module_mie \
  mie_sig.pyf \
  mie_source/anbn.f mie_source/devel.f mie_source/fichid.f \
  mie_source/gauleg.f mie_source/pitau.f mie_source/rminmax.f \
  mie_source/scatmat.f mie_source/sizedis.f mie_source/writsc.f

compile_module module_readmie \
  readmie_sig.pyf \
  readmie_source/read_mie_output.f readmie_source/file2coefs.f

compile_module module_mieshell \
  mieshell_sig.pyf \
  mieshell_source/bhcoat.f mieshell_source/gauleg.f mieshell_source/pitau.f \
  mieshell_source/scatmat.f mieshell_source/sizedis.f

compile_module module_dap \
  dap_sig.pyf \
  dap_source/adding.f dap_source/addlay.f dap_source/addsm.f \
  dap_source/assign.f dap_source/bmolecules.f dap_source/brack.f \
  dap_source/bstart.f dap_source/double.f dap_source/expbmu.f \
  dap_source/fillup.f dap_source/gauleg.f dap_source/init.f \
  dap_source/layer0.f dap_source/layerm.f dap_source/ldiapr.f \
  dap_source/newfou.f dap_source/nobot.f dap_source/notop.f \
  dap_source/ord1m.f dap_source/ord2m.f dap_source/prod.f \
  dap_source/rdiapr.f dap_source/renorm.f dap_source/scalzm.f \
  dap_source/setfou.f dap_source/setmu.f dap_source/setzm.f \
  dap_source/star.f dap_source/top2bot.f dap_source/trace.f \
  dap_source/transf.f dap_source/tstar.f

compile_module module_geos \
  geos_sig.pyf \
  geos_source/bracks.f geos_source/getgeos.f geos_source/rdfous.f \
  geos_source/read_dap.f geos_source/spline.f geos_source/splint.f \
  geos_source/wrout.f

# ── Smoke test ────────────────────────────────────────────────────────────────
echo ""
echo "── Smoke test ──────────────────────────────────────────"
python -c "
import module_mie, module_dap, module_geos, module_readmie, module_mieshell
import pymiedap.pymiedap as pmd
import pymiedap.ckdistribution
print('All modules imported successfully.')
m = pmd.Model()
print(f'pmd.Model() OK — default layers: {list(vars(m.layers).keys())}')
"

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "========================================================"
echo " Installation complete!"
echo ""
echo " To activate the environment in future sessions:"
echo "   source $VENV_ACTIVATE"
echo ""
echo " To run the example figures:"
echo "   cd $REPO_ROOT"
echo "   python examples/trees_stam_fig7.py --step dap"
echo "   python examples/trees_stam_fig7.py --step images"
echo "   python examples/trees_stam_fig7.py --step plot"
echo "========================================================"
