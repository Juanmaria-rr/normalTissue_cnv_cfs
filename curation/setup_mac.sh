#!/usr/bin/env bash
# One-time setup for the curation app on macOS / Linux.
# Creates a virtualenv inside this package and installs the dependencies.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$HERE/.venv"

echo "Curation app - setup"
echo "===================="
echo

# --- find a Python >= 3.11 -------------------------------------------------
# The app needs 3.11+; macOS ships 3.9, so look for a newer one explicitly
# before falling back to whatever `python3` happens to be.
PY=""
for cand in python3.13 python3.12 python3.11 python3; do
    if command -v "$cand" >/dev/null 2>&1; then
        if "$cand" -c 'import sys; sys.exit(0 if sys.version_info >= (3,11) else 1)' 2>/dev/null; then
            PY="$cand"
            break
        fi
    fi
done

if [ -z "$PY" ]; then
    echo "ERROR: no Python 3.11 or newer found." >&2
    echo >&2
    echo "macOS ships Python 3.9, which this app cannot use. Install a newer one:" >&2
    echo "    brew install python@3.11" >&2
    echo >&2
    echo "then run ./setup_mac.sh again." >&2
    if command -v python3 >/dev/null 2>&1; then
        echo >&2
        echo "(found python3 = $(python3 -V 2>&1))" >&2
    fi
    exit 1
fi

echo "Using $PY ($("$PY" -V 2>&1))"
echo

if [ -d "$VENV" ]; then
    echo "Virtualenv already exists at .venv - reusing it."
else
    echo "Creating virtualenv in .venv ..."
    "$PY" -m venv "$VENV"
fi

echo "Installing dependencies (this takes a minute) ..."
"$VENV/bin/python" -m pip install --quiet --upgrade pip
# Installed directly rather than via `pip install .[curation]`: installing the
# package would trigger setup.py, which compiles a C extension and needs Xcode
# command line tools. The curation app never uses that extension.
"$VENV/bin/python" -m pip install --quiet \
    "shiny>=1.0" "h5py>=3.10" "pandas>=2.1" "numpy>=1.26" \
    "scipy>=1.12" "scikit-learn>=1.3" "matplotlib>=3.8"

echo
echo "Verifying ..."
# $HERE goes in as argv[1] rather than being interpolated: the heredoc is
# quoted, so this works no matter which directory the script is invoked from.
"$VENV/bin/python" - "$HERE" <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1] + "/app")
import shiny, h5py, pandas, numpy, scipy, sklearn, matplotlib
from pylib.plotting.cn_profiles import plot_cn_profile
from pylib.plotting.grid_search import plot_grid_search_heatmap
from pylib.copynumber.noise import noise_estimate
print(f"  shiny {shiny.__version__}, h5py {h5py.__version__}, pandas {pandas.__version__}")
print("  curation library imports OK")
PYEOF

echo
echo "Setup complete. Launch the app with:"
echo "    ./run_curation.sh"
