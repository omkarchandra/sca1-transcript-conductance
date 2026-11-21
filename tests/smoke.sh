#!/usr/bin/env bash
# Parse and compile release sources without running analyses.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd -- "$ROOT"
PYTHON=${PYTHON:-python3}
RSCRIPT=${RSCRIPT:-Rscript}
# Keep bytecode outside the release tree and its checksum inventory.
export PYTHONPYCACHEPREFIX
PYTHONPYCACHEPREFIX=$(mktemp -d "${TMPDIR:-/tmp}/release-smoke-bytecode.XXXXXX")
while IFS= read -r -d '' file; do
    "$PYTHON" -m py_compile "$file"
done < <(find . -type f -name '*.py' -print0)
while IFS= read -r -d '' file; do
    bash -n "$file"
done < <(find . -type f \( -name '*.sh' \) -print0)
while IFS= read -r -d '' file; do
    "$RSCRIPT" --vanilla -e 'invisible(parse(file=commandArgs(TRUE)[1]))' "$file"
done < <(find . -type f -name '*.R' -print0)
"$PYTHON" config/paths.py --shell --config config/local.yaml.example > /dev/null
