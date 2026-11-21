#!/usr/bin/env bash
# Draws Figure 6 from the tables in inputs/.
# Usage (from this folder):  bash make_figure_6.sh OUT_DIR      OUT_DIR must not exist yet.
set -euo pipefail
OUT="${1:?usage: bash make_figure_6.sh OUT_DIR (a folder that does not exist yet)}"
mkdir "$OUT"
REF=../../data/reference_values/bushart2021_reference_values.md

python3 figure6_panels.py inputs/versions_without_synaptic_input.json "$REF" "$OUT/panels" \
    inputs/versions_with_synaptic_input.tsv

python3 figure6_page.py "$OUT/panels.png" "$OUT/panels.pdf" "$OUT/Figure_6"
