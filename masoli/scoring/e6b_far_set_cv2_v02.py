"""E6b v02: score the 37 °C far-set population.
v02 (never-run v01 f48cacf9 + provenance): skips the provenance/ folder, checks the kept list (b66e68e7) and each
run's factor_file_sha256 against provenance/<condition>.json.

Scorer v08 (bf8b41a3) is NOT used: its parse_condition (KI_RE, line 98) raises on every far_* condition name, and
Per run, from SPIKES.csv over [30000, 45000] ms
(both ends inclusive, as scorer v08 and E1): rate = spikes / 15 s; CV = SD (ddof 0) / mean ISI and CV2 (both need >= 2
ISIs, as scorer v08); pause fraction and log-ISI mode as the 23 °C driver v02 isi_metrics (copied from E1 v02).
WT settings = the 14 kept members' WT runs of the PRIMARY root, which must reproduce PER_RUN rate_hz and isi_cv exactly
(hard stop). Every far-set folder's SHA256SUMS and the SHA256SUMS_far_set_all are checked first; for every pair the
noise seed must equal the PRIMARY WT run's seed (20261013 + member). Read-only; writes only a NEW real_v01/.

  .venv-packet00/bin/python3 -I -B results/post_talk_20261024/E6b_far_set_cv2/e6b_far_set_cv2_v02.py --far-all-sha256 HASH
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import numpy as np

SSD = Path('${RESULTS_ROOT}')
SCORED = SSD / 'scored/masoli_pop_PRIMARY_45s_sig0.0065_v02_all840'
PRIMARY = SSD / 'masoli_pop_PRIMARY_45s_sig0.0065_v02'
FAR = SSD / 'masoli_pop_far_set_45s_sig0.0065_20261019'
FAR_ALL = 'SHA256SUMS_far_set_all'
KEPT_LIST = ('provenance/masoli_pop_far_set_kept_members_20261019.txt', 'b66e68e77158788a18b2c9c788e36523709f38feab303c3adbbb765f503ba7d5')
OUT = Path(__file__).resolve().parent / 'real_v01'
PINS = {SCORED / 'PER_RUN.tsv': 'b8b9f2a8d3f4e6e6a818626dab06f4bbbb6c2d0910f022b83fae1d7adf57ee5a',
        SCORED / 'PER_MEMBER.tsv': '148890d00ad1b45d2e5a04842d1083898c56242eedac1817b953ec3ba140d7fa',
        PRIMARY / 'SHA256SUMS_PRIMARY_all.txt': '37212b9ef7fc45d5584c82329c6514bd2fc601e0461f0ba959b1b7a8dc7bcedc'}
FAR_CONDITIONS = ('far_all_KI_nuclei_18w', 'far_wt_p90_18w', 'far_wt_p90_24w', 'far_wt_p90_30w',
                  'far_top_third_18w', 'far_top_third_24w', 'far_top_third_30w')
PRIMARY_REFERENCE = ('KI_18w_mid', 'KI_24w_mid', 'KI_30w_mid')   # all-KI-nuclei mid factors, beside the far sets
T0_MS, TSTOP_MS, WINDOW_S = 30000.0, 45000.0, 15.0
SEED_BASE = 20261013
PAUSE_FACTOR, MODE_BIN = 3.0, 0.05                               # 23 °C driver v02 lines 66-67
BANDS = {'cv': {'WT': (0.179, 0.350), 'SCA1': (0.314, 1.95)},    # driver v02 cv_percent bands / 100
         'cv2': {'WT': (0.176, 0.317), 'SCA1': (0.244, 0.614)}}  # driver v02 cv2 bands
NOTE = ('E6b exploratory scoring of the 37 °C far-set population (no new simulation); scorer v08 not used (rejects far_* '
        'names); E1 v02 metrics; window 30-45 s; WT = PRIMARY WT runs of the 14 kept sets.')


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_table(path, delimiter='\t'):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def check_sums(folder, sums):
    for line in Path(sums).read_text().splitlines():
        if line.strip():
            digest, name = line.split('  ', 1)
            if sha256(folder / name.removeprefix('./')) != digest:
                raise ValueError(f'checksum differs: {folder / name}')


def metrics(times):
    win = times[(times >= T0_MS) & (times <= TSTOP_MS)]
    isi = np.diff(win)
    out = dict(n_spikes=int(win.size), rate_hz=float(win.size / WINDOW_S), cv=float('nan'), cv2=float('nan'),
               pause_frac=float('nan'), log_isi_mode_ms=float('nan'))
    if isi.size >= 1:
        med = float(np.median(isi))
        out['pause_frac'] = float((isi > PAUSE_FACTOR * med).mean())
        edges = np.arange(-1.0, 5.0 + MODE_BIN / 2, MODE_BIN)
        hist, _ = np.histogram(np.log10(np.clip(isi, 10 ** edges[0], 10 ** edges[-1])), bins=edges)
        k = int(np.argmax(hist))
        out['log_isi_mode_ms'] = float(10 ** (0.5 * (edges[k] + edges[k + 1])))
    if isi.size >= 2:
        out['cv'] = float(np.std(isi, ddof=0) / np.mean(isi))
        out['cv2'] = float(np.mean(2.0 * np.abs(isi[1:] - isi[:-1]) / (isi[1:] + isi[:-1])))
    return out


def band_class(x, bands):
    if not np.isfinite(x):
        return 'undefined'
    wt, sca = bands['WT'][0] <= x <= bands['WT'][1], bands['SCA1'][0] <= x <= bands['SCA1'][1]
    return 'both_bands' if wt and sca else 'WT_band' if wt else 'SCA1_band' if sca else 'neither'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--far-all-sha256', required=True)
    args = ap.parse_args()
    for path, want in PINS.items():
        if sha256(path) != want:
            raise ValueError(f'input hash differs: {path}')
    if sha256(FAR / FAR_ALL) != args.far_all_sha256:
        raise ValueError("The SHA256SUMS_far_set_all differs from the reported hash")
    check_sums(FAR, FAR / FAR_ALL)

    per_run = read_table(SCORED / 'PER_RUN.tsv')
    kept = sorted({int(r['member_id']) for r in read_table(SCORED / 'PER_MEMBER.tsv') if r['kept'] == 'True'})
    if len(kept) != 14:
        raise ValueError('expected 14 kept sets')
    if sha256(FAR / KEPT_LIST[0]) != KEPT_LIST[1] or sorted(int(x) for x in (FAR / KEPT_LIST[0]).read_text().split()) != kept:
        raise ValueError("The kept-member list differs from PER_MEMBER's 14 kept sets")
    prim = {(int(r['member_id']), r['condition']): r for r in per_run}

    rows, wt_dev = [], 0.0
    for m in kept:
        for c in ('WT',) + PRIMARY_REFERENCE:
            r = prim[(m, c)]
            spikes = PRIMARY / r['folder'] / 'SPIKES.csv'
            if sha256(spikes) != r['spikes_sha256'] or int(r['noise_seed']) != SEED_BASE + m:
                raise ValueError(f'PRIMARY run differs from PER_RUN: {r["folder"]}')
            x = metrics(np.array([float(s['time_ms']) for s in read_table(spikes, ',')]))
            if c == 'WT':
                for key, col in (('rate_hz', 'rate_hz'), ('cv', 'isi_cv')):
                    old = float(r[col]) if r[col] not in ('', 'NA') else float('nan')
                    if not np.isfinite(old) or abs(old - x[key]) > 1e-12:
                        raise ValueError(f'E1 path disagrees with PRIMARY {col} on WT run {r["folder"]}: STOP')
                    wt_dev = max(wt_dev, abs(old - x[key]))
            rows.append(dict(source='PRIMARY', member_id=m, condition=c, folder=r['folder'], seed=int(r['noise_seed']), **x))
    folders = sorted(p for p in FAR.iterdir() if p.is_dir() and p.name != 'provenance')
    if len(folders) != 98:
        raise ValueError(f'expected 98 far-set folders, found {len(folders)}')
    seen = set()
    for f in folders:
        check_sums(f, f / 'SHA256SUMS')
        rec = json.loads((f / 'RECORD.json').read_text())
        m, c, seed = int(rec['member_id']), rec['condition'], int(rec['noise']['seed'])
        if m not in kept or c not in FAR_CONDITIONS or (m, c) in seen:
            raise ValueError(f'unexpected far-set folder: {f.name} {m} {c}')
        if seed != SEED_BASE + m or seed != int(prim[(m, 'WT')]['noise_seed']):
            raise ValueError(f'noise seed does not pair with the PRIMARY WT run: {f.name}')
        if rec['factor_file_sha256'] != sha256(FAR / 'provenance' / f'{c}.json'):
            raise ValueError(f'factor file differs from provenance copy: {f.name}')
        if float(rec['tstop_ms']) != TSTOP_MS:
            raise ValueError(f'tstop differs: {f.name}')
        seen.add((m, c))
        x = metrics(np.array([float(s['time_ms']) for s in read_table(f / 'SPIKES.csv', ',')]))
        rows.append(dict(source='far_set', member_id=m, condition=c, folder=f.name, seed=seed, **x))
    if len(seen) != 98:
        raise ValueError('far-set member x condition grid incomplete')

    by = {(r['member_id'], r['condition']): r for r in rows}
    summary = []
    for c in ('WT',) + PRIMARY_REFERENCE + FAR_CONDITIONS:
        rr = [by[(m, c)] for m in kept]
        s = dict(condition=c, n=len(rr), median_rate_hz=float(np.median([r['rate_hz'] for r in rr])))
        if c != 'WT':
            s['median_rate_ratio'] = float(np.median([r['rate_hz'] / by[(r['member_id'], 'WT')]['rate_hz'] for r in rr]))
        for metric in ('cv', 'cv2', 'pause_frac', 'log_isi_mode_ms'):
            v = np.array([r[metric] for r in rr])
            s[f'n_{metric}_defined'] = int(np.isfinite(v).sum())
            s[f'median_{metric}'] = float(np.nanmedian(v)) if np.isfinite(v).any() else float('nan')
            if c != 'WT':
                d = v - np.array([by[(r['member_id'], 'WT')][metric] for r in rr])
                s[f'median_change_{metric}'] = float(np.nanmedian(d)) if np.isfinite(d).any() else float('nan')
                s[f'n_{metric}_above_WT'] = int((d > 0).sum())
        for metric in ('cv', 'cv2'):
            classes = [band_class(r[metric], BANDS[metric]) for r in rr]
            for k in ('WT_band', 'SCA1_band', 'both_bands', 'neither', 'undefined'):
                s[f'{metric}_{k}'] = classes.count(k)
        summary.append(s)

    OUT.mkdir(parents=False, exist_ok=False)
    for name, table in (('E6b_PER_RUN.tsv', rows), ('E6b_CONDITION_SUMMARY.tsv', summary)):
        cols = list(dict.fromkeys(k for r in table for k in r))
        with (OUT / name).open('x', encoding='utf-8', newline='') as stream:
            stream.write(f'# {NOTE}\n')
            w = csv.DictWriter(stream, fieldnames=cols, delimiter='\t', lineterminator='\n', restval='')
            w.writeheader()
            w.writerows(table)
    record = {'note': NOTE, 'script_sha256': sha256(Path(__file__)), 'inputs': {str(p): h for p, h in PINS.items()},
              'far_all_sha256': args.far_all_sha256, 'far_folders_checked': len(folders), 'kept_members': kept,
              'wt_rate_cv_max_abs_dev_vs_PER_RUN': wt_dev, 'bands': BANDS, 'window_ms': [T0_MS, TSTOP_MS],
              'scorer_v08_used': False}
    with (OUT / 'RECORD.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(record, indent=1) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with (OUT / 'SHA256SUMS').open('x', encoding='utf-8') as stream:
        for p in files:
            stream.write(f'{sha256(p)}  {p.name}\n')
    print(json.dumps({'wt_dev': wt_dev, **{p.name: sha256(p) for p in files}}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
