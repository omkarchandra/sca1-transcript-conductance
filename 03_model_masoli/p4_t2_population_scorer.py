"""Population scorer for Track 2 pre-registration v2 (irregularity), v08 (v07 + Amendment S5a quantisation floor).

python -I -B 03_model_masoli/p4_t2_population_scorer.py          --runs-root DIR --t0-ms T0 --out DIR [--calibration]
Definitions follow. Code and
synthetic fixtures only; this script runs no simulation and reads only SPIKES.csv and RECORD.json.
v03: member_id (int) and condition come from RECORD.json; folder names (m000_c00, ...) are not parsed.
     Condition grammar as the files write it; per-condition reporting with verdicts only where the spec gives one.
v04: KI_18w_mid_grp_<label> conditions scored with the A1 reversal rule; prediction 3 (split 1) from BK_CaP and
     Kv3_Kv4_CaT; split 2 (BK_CaT, Kv3_Kv4_CaP) statistics only; calibration table sorted by noise_sigma_nA.
v05: the flag is --calibration; noise_sigma_nA, noise_tau_ms and noise_seed are read from RECORD.json when present
     (null otherwise) and carried into PER_RUN.tsv and CALIBRATION.tsv.
v06: noise settings are read from RECORD["noise"] = {sigma_nA, tau_ms, seed, sha256, ...} as the population
     driver v06 (2b57f435...) writes them; the top-level noise_sigma_nA / noise_tau_ms / noise_seed keys are read
     only when "noise" is absent; both forms present and disagreeing -> NoiseFieldConflict; sigma null or missing
     in both forms -> NoiseSigmaMissing(path), never defaulted to 0; noise.sha256 (null when the driver ran without
     noise) is carried into every per-run row; --calibration groups runs by the nested sigma_nA. Applies to both
     modes. The driver's other top-level keys (published, tail_isi_cv, amendment, ...) are not read.
v07: Amendment S5: the fixed CV-increase epsilon is replaced by a per-cell
     margin. Per run with >= 3 ISIs, isi_cv_se = SE of the ISI CV from a moving-block bootstrap (block length 10
     ISIs, or n if n < 10; blocks drawn with replacement at uniform start positions, concatenated, truncated to n;
     2,000 resamples; CV with ddof=0; SE = std of the 2,000 CVs with ddof=1; np.random.default_rng(20261013)
     created fresh per run). Fewer than 3 ISIs -> null. Per comparison cv_margin = 2*sqrt(se_cond^2 + se_wt^2)
     and cv_increased = cv_increase > cv_margin; a null SE on either side gives a null margin and no increase.
v08: Amendment S5a: a quantisation floor. Per run with >= 2 window ISIs,
     isi_q = (record_dt_ms / sqrt(6)) / mean(window ISIs) (sqrt(6): both ends of an ISI are rounded to the record
     grid); fewer ISIs -> null. record_dt_ms comes from RECORD.json (driver v06 writes 0.1); it must be a finite
     number > 0; missing or invalid in score mode -> RecordDtMissing(path) before any output; in calibration mode
     it is read when present (isi_q null otherwise). Per comparison cv_floor = sqrt(q_c^2 + q_wt^2) and
     cv_margin = max(2*sqrt(se_c^2 + se_wt^2), cv_floor); any of the four inputs null -> margin null, not increased.

Input layout: one run folder per (member, condition) directly under --runs-root, each holding
  SPIKES.csv   columns time_ms, peak_mV (header row required; may hold zero spikes)
  RECORD.json  with member_id (int), condition (str), tstop_ms (number)
Condition grammar:
  WT
  KI_<age>w_<low|mid|high>              kind ki        (mid: predictions; low/high: bounds, no verdict)
  KI_5w_mid_Hcn1_side                   kind hcn1_side (reported, no verdict)
  KI_<age>w_mid_restore_<gene>          kind restore   (A1 reversal verdict per gene)
  KI_<age>w_mid_grp_<label>             kind grp       (A1 reversal rule per label: BK_CaP, Kv3_Kv4_CaT, BK_CaT, Kv3_Kv4_CaP)
  KI_<age>w_mid_<label>                 kind sensitivity, any other suffix (descriptive only, no verdict)
Prediction 3 (split 1) holds if grp_BK_CaP is sufficient (share reversed >= 2/3) and grp_Kv3_Kv4_CaT is
not_sufficient (< 1/3); split 2 (grp_BK_CaT vs grp_Kv3_Kv4_CaP) gets the same statistics and no verdict.
Calibration mode (--calibration): one row per run folder with noise_sigma_nA (required in RECORD.json there;
member_id/condition optional), rate_hz and isi_cv (plus cv2, pauses, class) in CALIBRATION.tsv sorted by sigma;
no acceptance and no predictions; CALIBRATION_SUMMARY.json holds the input hashes; SHA256SUMS.

Per run, over the window [T0_MS, tstop_ms] (both ends inclusive, so a spike exactly at T0 or tstop counts):
  n_spikes, rate_hz = n_spikes / window_s
  isi_cv  = std(ISI, ddof=0) / mean(ISI)                          (needs >= 2 ISIs, else null)
  cv2     = mean of 2|ISI(i+1) - ISI(i)| / (ISI(i+1) + ISI(i))     (needs >= 2 ISIs, else null)
  pauses  = ISIs >= PAUSE_MS (count and durations); lead_gap_ms (T0 to first spike) and trail_gap_ms
            (last spike to tstop) are reported and not counted as pauses
  run_class = firing if n_spikes >= 1 else silent                  (no V trace: no plateau detection)
Population acceptance (spec section 4), WT runs only: kept if 20 <= rate_hz <= 120 and isi_cv < 0.2.
Per kept member and non-WT condition, against the member's own WT run:
  cv_increase = isi_cv(cond) - isi_cv(WT); rate_ratio = rate_hz(cond) / rate_hz(WT)
  a silent run has no isi_cv: cv_increase null (counts as no increase), rate_ratio 0
  isi_cv_se = bootstrap SE of isi_cv (see v07 above; needs >= 3 ISIs, else null)
  isi_q     = (record_dt_ms / sqrt(6)) / mean(window ISIs)     (needs >= 2 ISIs and record_dt_ms, else null)
  cv_floor  = sqrt(isi_q(cond)^2 + isi_q(WT)^2)                 (null if either q is null)
  cv_margin = max(2 * sqrt(isi_cv_se(cond)^2 + isi_cv_se(WT)^2), cv_floor)   (null if any input is null)
  cv_increased = cv_increase is not null and cv_margin is not null and cv_increase > cv_margin
  reversed (restore and grp) = |cv_increase| <= 0.05 and |rate_ratio - 1| <= 0.25
Per condition, over the kept members that have that run: n, n_cv_increased, share_cv_increased,
  median_cv_increase (nulls dropped), median_rate_ratio; then
  ki mid, age != 5:  prediction_1; change_detected = share >= 2/3 and median_cv_increase >= 0.10 and
                     |median_rate_ratio - 1| <= 0.25; expected change
  ki mid, age == 5:  prediction_2; same test; expected no change; consistent_with_expectation flag
  ki low/high, hcn1_side, sensitivity: statistics only, verdict null
  restore, grp:      A1 verdict: sufficient if share_reversed >= 2/3, not_sufficient if < 1/3, else partial
Outputs (exclusive create, nothing deleted): PER_RUN.tsv, PER_MEMBER.tsv, SUMMARY.json, SHA256SUMS
(score mode); CALIBRATION.tsv, CALIBRATION_SUMMARY.json, SHA256SUMS (calibrate mode).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np

PAUSE_MS = 500.0
RATE_MIN_HZ, RATE_MAX_HZ, CV_MAX = 20.0, 120.0, 0.2
SHARE_HOLD, CV_INCREASE_MIN, RATE_TOL = 2.0 / 3.0, 0.10, 0.25
REVERSAL_CV_TOL, SHARE_NOT_SUFFICIENT = 0.05, 1.0 / 3.0
BOOT_BLOCK, BOOT_B, BOOT_SEED = 10, 2000, 20261013
CV_INCREASE_RULE = ('delta > max(2*sqrt(se_c^2+se_wt^2), sqrt(q_c^2+q_wt^2)), q = (record_dt/sqrt(6))/mean ISI; '
                    'moving-block bootstrap, block 10, B 2000, seed 20261013')
SPLIT_1 = ('BK_CaP', 'Kv3_Kv4_CaT')      # prediction 3: first sufficient, second not_sufficient
SPLIT_2 = ('BK_CaT', 'Kv3_Kv4_CaP')      # statistics only, no verdict
CALIBRATION_COLUMNS = ['noise_sigma_nA', 'noise_tau_ms', 'noise_seed', 'noise_sha256', 'folder', 'member_id', 'condition',
                       'spikes_sha256', 'record_sha256', 't0_ms', 'tstop_ms', 'window_s', 'n_spikes', 'rate_hz',
                       'isi_cv', 'isi_cv_se', 'isi_q', 'cv2', 'n_pauses', 'run_class']
KI_RE = re.compile(r'^KI_(\d+)w_(low|mid|high)(?:_(.+))?$')
NOISE_KEYS = ('noise_sigma_nA', 'noise_tau_ms', 'noise_seed', 'noise_sha256')
NESTED_NOISE = {'noise_sigma_nA': 'sigma_nA', 'noise_tau_ms': 'tau_ms', 'noise_seed': 'seed', 'noise_sha256': 'sha256'}
RUN_COLUMNS = ['member_id', 'condition', 'kind', 'age', 'level', 'gene', 'label', 'folder', 'spikes_sha256',
               'record_sha256', 'noise_sigma_nA', 'noise_tau_ms', 'noise_seed', 'noise_sha256', 't0_ms', 'tstop_ms', 'window_s',
               'n_spikes_total', 'n_spikes', 'rate_hz', 'isi_cv', 'isi_cv_se', 'isi_q', 'cv2', 'n_pauses', 'pause_durations_ms',
               'lead_gap_ms', 'trail_gap_ms', 'run_class']
MEMBER_COLUMNS = ['member_id', 'kept', 'acceptance_reason', 'wt_rate_hz', 'wt_isi_cv', 'condition', 'kind',
                  'cv_increase', 'cv_floor', 'cv_margin', 'rate_ratio', 'cv_increased', 'reversed']


class NoiseFieldConflict(ValueError):
    """RECORD.json carries both the nested "noise" object and top-level noise_* keys, and they disagree."""


class RecordDtMissing(ValueError):
    """record_dt_ms is missing or not a finite number > 0 in RECORD.json (required in score mode)."""

    def __init__(self, path):
        super().__init__(f'{path}: record_dt_ms is missing or invalid in RECORD.json (needs a finite number > 0)')
        self.path = str(path)


class NoiseSigmaMissing(ValueError):
    """noise sigma is null or absent in both the nested and the top-level form; it is never defaulted to 0."""

    def __init__(self, path):
        super().__init__(f'{path}: noise sigma_nA is missing or null in RECORD.json (nested "noise" and top-level forms)')
        self.path = str(path)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_noise(record, path):
    """Return {noise_sigma_nA, noise_tau_ms, noise_seed, noise_sha256}; see the module docstring (v06)."""
    nested = record.get('noise')
    top = {k: record[k] for k in ('noise_sigma_nA', 'noise_tau_ms', 'noise_seed') if k in record}
    if nested is not None:
        if not isinstance(nested, dict):
            raise ValueError(f'{path}: RECORD.json "noise" must be an object')
        noise = {k: nested.get(n) for k, n in NESTED_NOISE.items()}
        for k, v in top.items():
            if v != noise[k]:
                raise NoiseFieldConflict(f'{path}: {k}={v!r} (top level) disagrees with noise.{NESTED_NOISE[k]}='
                                         f'{noise[k]!r} (nested)')
    else:
        noise = {k: top.get(k) for k in ('noise_sigma_nA', 'noise_tau_ms', 'noise_seed')}
        noise['noise_sha256'] = None
    if noise['noise_sigma_nA'] is None:
        raise NoiseSigmaMissing(path)
    for key in ('noise_sigma_nA', 'noise_tau_ms'):
        if noise[key] is not None:
            if isinstance(noise[key], bool):
                raise ValueError(f'{path}: {key} must be a number')
            noise[key] = float(noise[key])
            if not np.isfinite(noise[key]):
                raise ValueError(f'{path}: {key} must be finite')
    seed = noise['noise_seed']
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
        raise ValueError(f'{path}: noise seed must be an int')
    sha = noise['noise_sha256']
    if sha is not None and not (isinstance(sha, str) and re.fullmatch(r'[0-9a-f]{64}', sha)):
        raise ValueError(f'{path}: noise sha256 must be a 64-hex string or null')
    return noise


def fmt(x):
    if x is None:
        return 'NA'
    if isinstance(x, float):
        return repr(x)
    if isinstance(x, list):
        return ';'.join(fmt(v) for v in x)
    return str(x)


def parse_condition(condition: str):
    """Return dict(kind, age, level, gene, label)."""
    if condition == 'WT':
        return dict(kind='wt', age=None, level=None, gene=None, label=None)
    m = KI_RE.match(condition)
    if not m:
        raise ValueError(f'condition {condition!r} does not match WT | KI_<age>w_<low|mid|high>[_suffix]')
    age, level, suffix = int(m.group(1)), m.group(2), m.group(3)
    if suffix is None:
        return dict(kind='ki', age=age, level=level, gene=None, label=None)
    if suffix == 'Hcn1_side':
        return dict(kind='hcn1_side', age=age, level=level, gene='Hcn1', label=suffix)
    if suffix.startswith('restore_') and len(suffix) > len('restore_'):
        return dict(kind='restore', age=age, level=level, gene=suffix[len('restore_'):], label=suffix)
    if suffix.startswith('grp_') and len(suffix) > len('grp_'):
        return dict(kind='grp', age=age, level=level, gene=None, label=suffix[len('grp_'):])
    return dict(kind='sensitivity', age=age, level=level, gene=None, label=suffix)


def read_run(folder: Path, calibrate: bool = False):
    spikes_bytes = (folder / 'SPIKES.csv').read_bytes()
    record_bytes = (folder / 'RECORD.json').read_bytes()
    record = json.loads(record_bytes)
    required = ('tstop_ms',) if calibrate else ('member_id', 'condition', 'tstop_ms')
    for key in required:
        if key not in record:
            raise ValueError(f'{folder.name}: RECORD.json lacks {key}')
    member_id = record.get('member_id')
    if member_id is not None and (isinstance(member_id, bool) or not isinstance(member_id, int)):
        raise ValueError(f'{folder.name}: member_id must be an int')
    condition = None if record.get('condition') is None else str(record['condition'])
    noise = read_noise(record, folder)
    record_dt = record.get('record_dt_ms')
    if record_dt is None and calibrate:
        pass                                   # optional in calibration mode: isi_q stays null
    else:
        if (record_dt is None or isinstance(record_dt, bool) or not isinstance(record_dt, (int, float))
                or not np.isfinite(float(record_dt)) or float(record_dt) <= 0):
            raise RecordDtMissing(folder)
        record_dt = float(record_dt)
    tstop = float(record['tstop_ms'])
    if not np.isfinite(tstop) or tstop <= 0:
        raise ValueError(f'{folder.name}: tstop_ms must be positive and finite')
    rows = list(csv.DictReader(spikes_bytes.decode().splitlines()))
    if rows and not {'time_ms', 'peak_mV'} <= set(rows[0]):
        raise ValueError(f'{folder.name}: SPIKES.csv needs columns time_ms, peak_mV')
    times = np.array([float(r['time_ms']) for r in rows], dtype=float)
    if times.size and (not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0)):
        raise ValueError(f'{folder.name}: spike times must be finite and strictly increasing')
    if times.size and times[-1] > tstop:
        raise ValueError(f'{folder.name}: spike after tstop_ms')
    parsed = (parse_condition(condition) if condition is not None and not calibrate else
              dict(kind=None, age=None, level=None, gene=None, label=None))
    return dict(member_id=member_id, condition=condition, parsed=parsed, folder=folder.name, noise=noise,
                record_dt_ms=record_dt, tstop_ms=tstop, times=times, spikes_sha256=digest(spikes_bytes),
                record_sha256=digest(record_bytes))


def bootstrap_cv_se(isi):
    """Moving-block bootstrap SE of the ISI CV (Amendment S5). Needs >= 3 ISIs; see the module docstring."""
    n = int(isi.size)
    if n < 3:
        return None
    block = min(BOOT_BLOCK, n)
    n_blocks = -(-n // block)
    rng = np.random.default_rng(BOOT_SEED)
    starts = rng.integers(0, n - block + 1, size=(BOOT_B, n_blocks))
    idx = (starts[:, :, None] + np.arange(block)[None, None, :]).reshape(BOOT_B, n_blocks * block)[:, :n]
    samples = isi[idx]
    cvs = np.std(samples, axis=1, ddof=0) / np.mean(samples, axis=1)
    return float(np.std(cvs, ddof=1))


def score_run(run, t0_ms: float):
    tstop = run['tstop_ms']
    if t0_ms < 0 or t0_ms >= tstop:
        raise ValueError(f"{run['folder']}: T0_MS {t0_ms} must satisfy 0 <= T0 < tstop_ms {tstop}")
    window_s = (tstop - t0_ms) / 1000.0
    times = run['times']
    inside = times[(times >= t0_ms) & (times <= tstop)]
    n = int(inside.size)
    isi = np.diff(inside)
    cv = cv2 = None
    if isi.size >= 2:
        cv = float(np.std(isi, ddof=0) / np.mean(isi))
        cv2 = float(np.mean(2.0 * np.abs(isi[1:] - isi[:-1]) / (isi[1:] + isi[:-1])))
    cv_se = bootstrap_cv_se(isi)
    q = (None if isi.size < 2 or run['record_dt_ms'] is None else
         float((run['record_dt_ms'] / np.sqrt(6.0)) / np.mean(isi)))
    pauses = [float(d) for d in isi[isi >= PAUSE_MS]]
    p = run['parsed']
    return dict(member_id=run['member_id'], condition=run['condition'], kind=p['kind'], age=p['age'],
                level=p['level'], gene=p['gene'], label=p['label'], folder=run['folder'],
                **run['noise'],
                spikes_sha256=run['spikes_sha256'], record_sha256=run['record_sha256'],
                t0_ms=float(t0_ms), tstop_ms=float(tstop), window_s=float(window_s),
                n_spikes_total=int(times.size), n_spikes=n, rate_hz=float(n / window_s),
                isi_cv=cv, isi_cv_se=cv_se, isi_q=q, cv2=cv2, n_pauses=len(pauses), pause_durations_ms=pauses,
                lead_gap_ms=float(inside[0] - t0_ms) if n else float(tstop - t0_ms),
                trail_gap_ms=float(tstop - inside[-1]) if n else float(tstop - t0_ms),
                run_class='firing' if n >= 1 else 'silent')


def is_kept(wt):
    if wt is None:
        return False, 'no_WT_run'
    if wt['isi_cv'] is None:
        return False, 'WT_fewer_than_3_spikes'
    if not RATE_MIN_HZ <= wt['rate_hz'] <= RATE_MAX_HZ:
        return False, 'WT_rate_outside_20_120_Hz'
    if not wt['isi_cv'] < CV_MAX:
        return False, 'WT_isi_cv_not_below_0.2'
    return True, 'kept'


def median_or_none(values):
    values = [v for v in values if v is not None]
    return float(np.median(values)) if values else None


def summarise(scored):
    by_member = {}
    for r in scored:
        conds = by_member.setdefault(r['member_id'], {})
        if r['condition'] in conds:
            raise ValueError(f"duplicate run for member {r['member_id']} condition {r['condition']}")
        conds[r['condition']] = r
    members = []
    for member_id in sorted(by_member):
        conds = by_member[member_id]
        wt = conds.get('WT')
        kept, reason = is_kept(wt)
        entry = dict(member_id=member_id, kept=kept, acceptance_reason=reason,
                     wt_rate_hz=None if wt is None else wt['rate_hz'],
                     wt_isi_cv=None if wt is None else wt['isi_cv'], comparisons={})
        if kept:
            for cond, r in sorted(conds.items()):
                if r['kind'] == 'wt':
                    continue
                ratio = float(r['rate_hz'] / wt['rate_hz'])
                cv_inc = None if r['isi_cv'] is None else float(r['isi_cv'] - wt['isi_cv'])
                floor = (None if r['isi_q'] is None or wt['isi_q'] is None else
                         float(np.sqrt(r['isi_q'] ** 2 + wt['isi_q'] ** 2)))
                margin = (None if r['isi_cv_se'] is None or wt['isi_cv_se'] is None or floor is None else
                          float(max(2.0 * np.sqrt(r['isi_cv_se'] ** 2 + wt['isi_cv_se'] ** 2), floor)))
                increased = bool(cv_inc is not None and margin is not None and cv_inc > margin)
                reversed_ = None
                if r['kind'] in ('restore', 'grp'):
                    reversed_ = bool(cv_inc is not None and abs(cv_inc) <= REVERSAL_CV_TOL
                                     and abs(ratio - 1.0) <= RATE_TOL)
                entry['comparisons'][cond] = dict(kind=r['kind'], age=r['age'], level=r['level'], gene=r['gene'],
                                                  label=r['label'], cv_increase=cv_inc, cv_floor=floor, cv_margin=margin,
                                                  rate_ratio=ratio,
                                                  cv_increased=increased, reversed=reversed_,
                                                  run_class=r['run_class'])
        members.append(entry)
    kept_members = [m for m in members if m['kept']]
    conditions = {}
    for cond in sorted({c for m in kept_members for c in m['comparisons']}):
        rows = [m['comparisons'][cond] for m in kept_members if cond in m['comparisons']]
        p = parse_condition(cond)
        n_inc = int(sum(r['cv_increased'] for r in rows))
        share = float(n_inc / len(rows))
        med_cv = median_or_none([r['cv_increase'] for r in rows])
        med_ratio = median_or_none([r['rate_ratio'] for r in rows])
        out = dict(kind=p['kind'], age=p['age'], level=p['level'], gene=p['gene'], label=p['label'],
                   n_kept_members_with_run=len(rows), n_cv_increased=n_inc, share_cv_increased=share,
                   median_cv_increase=med_cv, median_rate_ratio=med_ratio,
                   n_silent=int(sum(r['run_class'] == 'silent' for r in rows)),
                   role=None, verdict=None)
        if p['kind'] == 'ki' and p['level'] == 'mid':
            holds = bool(share >= SHARE_HOLD and med_cv is not None and med_cv >= CV_INCREASE_MIN
                         and abs(med_ratio - 1.0) <= RATE_TOL)
            which = 'prediction_2' if p['age'] == 5 else 'prediction_1'
            out.update(role=which, change_detected=holds,
                       expected='no_change' if which == 'prediction_2' else 'change',
                       consistent_with_expectation=bool(holds == (which == 'prediction_1')),
                       verdict='change_detected' if holds else 'no_change_detected')
        elif p['kind'] == 'ki':
            out.update(role=f"bounds_{p['level']}")
        elif p['kind'] == 'hcn1_side':
            out.update(role='hcn1_side_report')
        elif p['kind'] in ('restore', 'grp'):
            n_rev = int(sum(bool(r['reversed']) for r in rows))
            share_rev = float(n_rev / len(rows))
            out.update(role='A1_reversal' if p['kind'] == 'restore' else 'grp_reversal',
                       n_reversed=n_rev, share_reversed=share_rev,
                       verdict=('sufficient' if share_rev >= SHARE_HOLD else
                                'not_sufficient' if share_rev < SHARE_NOT_SUFFICIENT else 'partial'))
        else:
            out.update(role='sensitivity_descriptive')
        conditions[cond] = out
    return members, conditions


def split_tests(conditions):
    """Prediction 3 (split 1) and the no-verdict split 2, from the grp_<label> condition verdicts."""
    grp = {v['label']: dict(condition=c, verdict=v['verdict'], n_reversed=v['n_reversed'],
                            n_kept_members_with_run=v['n_kept_members_with_run'], share_reversed=v['share_reversed'])
           for c, v in conditions.items() if v['role'] == 'grp_reversal'}
    a, b = SPLIT_1
    holds = None
    if a in grp and b in grp:
        holds = bool(grp[a]['verdict'] == 'sufficient' and grp[b]['verdict'] == 'not_sufficient')
    split_1 = dict(role='prediction_3', rule=f'{a} sufficient and {b} not_sufficient', labels=list(SPLIT_1),
                   groups={k: grp.get(k) for k in SPLIT_1}, holds=holds,
                   verdict=None if holds is None else ('holds' if holds else 'does_not_hold'),
                   note=None if holds is not None else 'a split-1 condition is missing among the kept members')
    split_2 = dict(role='split_2_descriptive', labels=list(SPLIT_2), groups={k: grp.get(k) for k in SPLIT_2},
                   holds=None, verdict=None, note='statistics only; no verdict by spec')
    return dict(split_1=split_1, split_2=split_2)


def member_rows(members):
    rows = []
    for m in members:
        base = dict(member_id=m['member_id'], kept=m['kept'], acceptance_reason=m['acceptance_reason'],
                    wt_rate_hz=m['wt_rate_hz'], wt_isi_cv=m['wt_isi_cv'])
        if not m['comparisons']:
            rows.append(dict(base, condition='', kind='', cv_increase=None, cv_floor=None, cv_margin=None,
                             rate_ratio=None, cv_increased=None, reversed=None))
        for cond, r in m['comparisons'].items():
            rows.append(dict(base, condition=cond, kind=r['kind'], cv_increase=r['cv_increase'],
                             cv_floor=r['cv_floor'], cv_margin=r['cv_margin'], rate_ratio=r['rate_ratio'],
                             cv_increased=r['cv_increased'], reversed=r['reversed']))
    return rows


def write_tsv(path: Path, columns, rows):
    with path.open('x') as f:
        f.write('\t'.join(columns) + '\n')
        for r in rows:
            f.write('\t'.join(fmt(r[c]) for c in columns) + '\n')


def calibrate(args, folders):
    runs = [read_run(p, calibrate=True) for p in folders]   # every input read and hashed before any output
    scored = [score_run(r, args.t0_ms) for r in runs]
    scored.sort(key=lambda r: (r['noise_sigma_nA'], r['folder']))
    args.out.mkdir(parents=False, exist_ok=False)
    write_tsv(args.out / 'CALIBRATION.tsv', CALIBRATION_COLUMNS, scored)
    groups = {}
    for r in scored:
        groups.setdefault(repr(r['noise_sigma_nA']), []).append(r)
    group_stats = {sigma: dict(noise_sigma_nA=rows[0]['noise_sigma_nA'], n_runs=len(rows),
                               folders=[r['folder'] for r in rows],
                               median_rate_hz=median_or_none([r['rate_hz'] for r in rows]),
                               median_isi_cv=median_or_none([r['isi_cv'] for r in rows]),
                               n_silent=int(sum(r['run_class'] == 'silent' for r in rows)))
                   for sigma, rows in groups.items()}
    summary = dict(mode='calibration', script=Path(__file__).name, script_sha256=digest(Path(__file__).read_bytes()),
                   runs_root=str(args.runs_root), t0_ms=float(args.t0_ms), n_runs=len(scored),
                   sorted_by='noise.sigma_nA ascending, then folder name', groups_by_sigma=group_stats,
                   rows=[{k: r[k] for k in CALIBRATION_COLUMNS} for r in scored],
                   input_sha256={r['folder']: dict(**{k: r[k] for k in NOISE_KEYS}, SPIKES_csv=r['spikes_sha256'],
                                                   RECORD_json=r['record_sha256']) for r in scored})
    with (args.out / 'CALIBRATION_SUMMARY.json').open('x') as f:
        f.write(json.dumps(summary, indent=2) + '\n')
    names = ['CALIBRATION.tsv', 'CALIBRATION_SUMMARY.json']
    with (args.out / 'SHA256SUMS').open('x') as f:
        f.write(''.join(f'{digest((args.out / n).read_bytes())}  {n}\n' for n in names))
    print(json.dumps(dict(out=str(args.out), mode='calibration', n_runs=len(scored),
                          sigma_order=[r['noise_sigma_nA'] for r in scored],
                          groups={k: v['n_runs'] for k, v in group_stats.items()}), indent=2))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--runs-root', required=True, type=Path)
    ap.add_argument('--t0-ms', required=True, type=float)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--calibration', action='store_true',
                    help='calibration table sorted by noise_sigma_nA; no acceptance, no predictions')
    args = ap.parse_args(argv)
    folders = sorted(p for p in args.runs_root.iterdir() if p.is_dir() and (p / 'SPIKES.csv').exists())
    if not folders:
        raise SystemExit(f'no run folders with SPIKES.csv under {args.runs_root}')
    if args.calibration:
        calibrate(args, folders)
        return
    runs = [read_run(p) for p in folders]          # every input read and hashed before any output
    scored = [score_run(r, args.t0_ms) for r in runs]
    scored.sort(key=lambda r: (r['member_id'], r['condition']))
    members, conditions = summarise(scored)
    splits = split_tests(conditions)
    args.out.mkdir(parents=False, exist_ok=False)   # exclusive: never writes into an existing folder
    write_tsv(args.out / 'PER_RUN.tsv', RUN_COLUMNS, scored)
    write_tsv(args.out / 'PER_MEMBER.tsv', MEMBER_COLUMNS, member_rows(members))
    summary = dict(
        script=Path(__file__).name, script_sha256=digest(Path(__file__).read_bytes()),
        runs_root=str(args.runs_root), t0_ms=float(args.t0_ms),
        parameters=dict(pause_ms=PAUSE_MS, rate_min_hz=RATE_MIN_HZ, rate_max_hz=RATE_MAX_HZ, cv_max=CV_MAX,
                        share_hold=SHARE_HOLD, cv_increase_min=CV_INCREASE_MIN, rate_tol=RATE_TOL,
                        reversal_cv_tol=REVERSAL_CV_TOL, share_not_sufficient=SHARE_NOT_SUFFICIENT,
                        cv_increase_rule=CV_INCREASE_RULE, bootstrap=dict(block=BOOT_BLOCK, B=BOOT_B, seed=BOOT_SEED)),
        n_runs=len(scored), n_members=len(members), n_kept=int(sum(m['kept'] for m in members)),
        kept_member_ids=[m['member_id'] for m in members if m['kept']],
        rejected_members=[dict(member_id=m['member_id'], reason=m['acceptance_reason'])
                          for m in members if not m['kept']],
        run_classes={r['folder']: r['run_class'] for r in scored},
        predictions={c: v for c, v in conditions.items() if v['role'] in ('prediction_1', 'prediction_2')},
        bounds={c: v for c, v in conditions.items() if v['role'] in ('bounds_low', 'bounds_high')},
        hcn1_side={c: v for c, v in conditions.items() if v['role'] == 'hcn1_side_report'},
        reversal={c: v for c, v in conditions.items() if v['role'] == 'A1_reversal'},
        grp_reversal={c: v for c, v in conditions.items() if v['role'] == 'grp_reversal'},
        prediction_3=splits['split_1'], split_2=splits['split_2'],
        sensitivity={c: v for c, v in conditions.items() if v['role'] == 'sensitivity_descriptive'},
        members=members,
        input_sha256={r['folder']: dict(member_id=r['member_id'], condition=r['condition'],
                                        SPIKES_csv=r['spikes_sha256'], RECORD_json=r['record_sha256'],
                                        **{k: r[k] for k in NOISE_KEYS})
                      for r in scored})
    with (args.out / 'SUMMARY.json').open('x') as f:
        f.write(json.dumps(summary, indent=2) + '\n')
    names = ['PER_RUN.tsv', 'PER_MEMBER.tsv', 'SUMMARY.json']
    with (args.out / 'SHA256SUMS').open('x') as f:
        f.write(''.join(f'{digest((args.out / n).read_bytes())}  {n}\n' for n in names))
    print(json.dumps(dict(out=str(args.out), n_runs=len(scored), n_kept=summary['n_kept'],
                          verdicts={c: v['verdict'] for c, v in conditions.items() if v['verdict'] is not None},
                          prediction_3=splits['split_1']['verdict']),
                     indent=2))


if __name__ == '__main__':
    main()

#%% debug
