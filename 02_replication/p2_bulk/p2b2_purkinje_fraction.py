"""v02: the gate allows the 6-significant-digit rounding of P2b_TABLE.tsv (v01 used 1e-12 and stopped at 4.8e-7).
P2b-2 (descriptive): P2b with Kcnmb4 and Cacna2d2 added..

Takes P2b's table (f, observed) unchanged and fills the 12-week Purkinje log2FC of all nine genes from E3's LIMMA_RESULTS (now on
the SSD); gate: the seven genes P2b took from E3's TWELVE_WEEK_GENES must agree with LIMMA_RESULTS to 1e-12. Predicted bulk =
f x log2FC; Spearman rho over the five (unchanged) and over all nine. Read-only; writes a NEW real_v01/ next to this script.

  python -I -B 02_replication/p2_bulk/p2b2_purkinje_fraction.py
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/p2b2_purkinje_fraction'))
P2B = ROOT / _release_paths.expand('${DATA_ROOT}/publication/P2b_purkinje_fraction/real_v01/P2b_TABLE.tsv')
LIMMA = Path(_release_paths.expand('${REPO_ROOT}/data/snrnaseq/E3_LIMMA_RESULTS.tsv'))
PINS = {P2B: 'd3c2d34ba1d1750c0997a14e3e71a5678f3eccfb2eb7c82a3ce3542cb33de45a',
        LIMMA: 'eab5ccbedff4c8b340ed6b6a9c43920bcfe7aa308387a2702f3cd8ab4c101e67'}
FIVE = ('Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g')
NOTE = 'P2b-2 (descriptive): P2b with all nine 12-week Purkinje log2FC from E3 LIMMA_RESULTS; f and observed from P2b unchanged.'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def main():
    for p, want in PINS.items():
        if sha(p) != want:
            raise ValueError(f'input hash differs: {p}')
    t = pd.read_csv(P2B, sep='\t', comment='#')
    lim = pd.read_csv(LIMMA, sep='\t')
    lim = lim[(lim.population == 'Purkinje cell') & (lim.age_weeks == 12)].set_index('dataset_symbol')
    e3 = lim.loc[t.gene]
    had = t.purkinje_12w_log2fc.notna()
    old, new = t.purkinje_12w_log2fc[had].to_numpy(), e3.log2FC[had.to_numpy()].to_numpy()
    dev = float(np.abs(old - new).max())
    if np.any(np.abs(old - new) > 5e-6 * np.abs(new)):        # v02: P2b_TABLE stores 6 significant digits (%.6g)
        raise ValueError(f'gate: P2b log2FC differs from LIMMA_RESULTS by {dev}')
    t['purkinje_12w_log2fc'] = e3.log2FC.to_numpy()
    t['purkinje_12w_CI_L'], t['purkinje_12w_CI_R'], t['purkinje_12w_p'] = e3.CI_L.to_numpy(), e3.CI_R.to_numpy(), e3.p.to_numpy()
    t['predicted_bulk_log2fc'] = t.f_purkinje * t.purkinje_12w_log2fc
    t['note'] = t.note.fillna('').str.replace('not in E3 TWELVE_WEEK_GENES; E3 LIMMA_RESULTS not on local disk; left out', 'log2FC from E3 LIMMA_RESULTS (P2b-2)')
    r5 = spearmanr(t[t.gene.isin(FIVE)].predicted_bulk_log2fc, t[t.gene.isin(FIVE)].observed_log2fc).statistic
    r9 = spearmanr(t.predicted_bulk_log2fc, t.observed_log2fc).statistic
    OUT.mkdir(parents=False, exist_ok=False)
    buf = io.StringIO(); buf.write('# ' + NOTE + '\n'); t.to_csv(buf, sep='\t', index=False, float_format='%.6g')
    with open(OUT / 'P2b2_TABLE.tsv', 'x') as f:
        f.write(buf.getvalue())
    rec = dict(note=NOTE, task='P2b-2', script_sha256=sha(Path(__file__)), inputs={str(p): h for p, h in PINS.items()},
               gate_max_abs_vs_P2b=dev, rho_five=float(r5), rho_nine=float(r9), finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print(t[['gene', 'f_purkinje', 'purkinje_12w_log2fc', 'predicted_bulk_log2fc', 'observed_log2fc']].round(3).to_string(index=False))
    print(json.dumps(dict(rho_five=r5, rho_nine=r9, gate=dev, **{p.name: sha(p) for p in files}), indent=1, default=float))


if __name__ == '__main__':
    main()

#%% debug
