"""Purkinje-only counts for the manuscript spread paragraph (

Inputs (pinned): the accepted A1 within-mouse spread tables (TREND_BY_GENOTYPE_PC.tsv cbbebf8e, KI_VS_WT_BY_AGE_PC.tsv a443708b) and the accepted
A3 bimodality calls (LANE_AGE_CALLS.tsv 2a2a9bd6). Slopes as in the P4 reading (section 3): a straight-line fit per week to the genotype-mean SD;
primary = the 4 discovery ages (5/18/24/30 W), and the A1 table's saved 5-age slope is reported alongside. The same counts over all 12 components
are recomputed as a check against the paragraph's earlier numbers (10, 9, 12, 11 and 7 of 12). Writes a NEW real_v01/SPREAD_COUNTS.json + SHA256SUMS.

  python -I -B 01_snrnaseq/spread/spread_counts_purkinje.py
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
import json

import numpy as np
import pandas as pd

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/spread_counts_purkinje'))
TREND = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_lane_tightening_A1/TREND_BY_GENOTYPE_PC.tsv')
BYAGE = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_lane_tightening_A1/KI_VS_WT_BY_AGE_PC.tsv')
BIMOD = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_a3_bimodality/LANE_AGE_CALLS.tsv')
PINS = {TREND: 'cbbebf8e8dd8a2b59730da3c8822ee0886192cf796920dab5b06ec9bda0cf087', BYAGE: 'a443708bf8eeba10a9d55a8a3eb07c00820da9cf3a1b6c61110e293c3dbbbdde',
        BIMOD: '2a2a9bd6d3ca0be8c149bcf44ab6db8996846e41b058c9100d426f03df280880'}
DISC = (5, 18, 24, 30); ALL5 = (5, 12, 18, 24, 30)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def counts(t, k30):
    s = {}
    for name, ages in (('discovery_4_ages_primary', DISC), ('saved_5_ages', ALL5)):
        sl = t.apply(lambda r: np.polyfit(ages, [r[f'mean_sd_{a}W'] for a in ages], 1)[0], axis=1)
        w = pd.DataFrame(dict(PC=t.PC, genotype=t.genotype, population=t.population, slope=sl)).pivot_table(index=['population', 'PC'], columns='genotype', values='slope')
        s[name] = dict(n_components=len(w), KI_spread_rising=int((w.KI > 0).sum()), WT_spread_falling=int((w.WT < 0).sum()), KI_slope_greater_than_WT=int((w.KI > w.WT).sum()))
    s['saved_5_ages']['matches_saved_slope_column'] = bool(np.allclose(t.apply(lambda r: np.polyfit(ALL5, [r[f'mean_sd_{a}W'] for a in ALL5], 1)[0], axis=1), t.slope_sd_per_week))
    s['age_30W'] = dict(KI_spread_greater_than_WT=int((k30.diff_KI_minus_WT > 0).sum()), of_which_complete_separation_KI_higher=int(((k30.diff_KI_minus_WT > 0) & (k30.separation == 'KI_higher')).sum()),
                        complete_separation_KI_lower=int((k30.separation == 'KI_lower').sum()))
    return s


def main():
    for p, want in PINS.items():
        if sha(p) != want:
            raise ValueError(f'input hash differs: {p}')
    t = pd.read_csv(TREND, sep='\t'); a = pd.read_csv(BYAGE, sep='\t'); b = pd.read_csv(BIMOD, sep='\t')
    k30 = a[a.age_weeks == 30]
    pk = t.population == 'Purkinje_cell'; assert sorted(t[pk].PC.unique()) == [5, 6, 11, 30]
    res = dict(purkinje=counts(t[pk].reset_index(drop=True), k30[k30.population == 'Purkinje_cell']), all_12_check=counts(t, k30))
    bp = b[b.population == 'Purkinje_cell']; assert len(bp) == 8
    res['purkinje']['bimodality_24_30W'] = dict(n_component_ages=len(bp), KI_component_ages_with_any_two_group_mouse=int((bp.KI_flagged > 0).sum()),
                                                WT_component_ages_with_any_two_group_mouse=int((bp.WT_flagged > 0).sum()), KI_specific_two_group_meeting_rule=int((bp.KI_specific_two_group == 'yes').sum()),
                                                flagged_rows=bp[(bp.KI_flagged > 0) | (bp.WT_flagged > 0)][['PC', 'age_weeks', 'KI_flagged', 'WT_flagged']].to_dict('records'))
    res.update(inputs={str(p): h for p, h in PINS.items()}, script_sha256=sha(Path(__file__)), finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    OUT.mkdir(parents=False, exist_ok=False)
    with open(OUT / 'SPREAD_COUNTS.json', 'x') as f:
        f.write(json.dumps(res, indent=1, default=int) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        f.write(f"{sha(OUT / 'SPREAD_COUNTS.json')}  SPREAD_COUNTS.json\n")
    print(json.dumps(res, indent=1, default=int))


if __name__ == '__main__':
    main()

#%% debug
