#!/bin/bash
# P14 environment. = p13_venv_v02.sh (55438547...) with BASE changed.
# Rebuilds the stage-2 environment: python 3.10.1 module, NEURON 8.2.6 from the cached wheel,
# numpy 1.26.4, find-libpython 0.5.1, packaging 26.3 (all as pip82_install.log of 26 Oct). v02 = v01 (b89d92ca...) with packaging pinned. Exclusive create; no retry.
set -euo pipefail
# Needs Python 3.10.1.
BASE=${PROJECT_ROOT}/p14_zang2021_20261110
VENV=$BASE/venv
test -d "$BASE/inputs" && test ! -e "$VENV"
python3 -m venv "$VENV"
"$VENV/bin/python" -m pip install --no-cache-dir numpy==1.26.4 find-libpython==0.5.1 packaging==26.3 > "$BASE/logs/venv_pip1.log" 2>&1
"$VENV/bin/python" -m pip install --no-cache-dir --no-deps \
    "$BASE/inputs/wheel/NEURON-8.2.6-cp310-cp310-manylinux_2_17_x86_64.manylinux2014_x86_64.whl" > "$BASE/logs/venv_pip2.log" 2>&1
"$VENV/bin/python" -I -B -c 'import neuron, numpy, sys; print(sys.version.split()[0], neuron.__version__, numpy.__version__)' \
    | tee "$BASE/receipts/VENV_VERSIONS.txt"
"$VENV/bin/python" -m pip freeze > "$BASE/receipts/VENV_FREEZE.txt"
grep -q '^3.10.1 8.2.6' "$BASE/receipts/VENV_VERSIONS.txt" && grep -q ' 1.26.4$' "$BASE/receipts/VENV_VERSIONS.txt"
echo VENV_OK
