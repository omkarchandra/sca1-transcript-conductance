""" independent re-score of the Track 2 Stage A traces (pre-specified).

python -I -B 05_model_forrest/stageA_rescore.py
Reads only V_TRACE.tsv, V_DEND_TRACE.tsv and SPIKES.tsv of the two pinned Stage A runs (hashes checked
before use), never the RECORD.json, SHA256SUMS or run reports. Scoring is written from review's spec, not
from the code. Writes only under ${DATA_ROOT}/analysis/stageA_rescore/ with exclusive create.
No simulation, nothing deleted.

Definitions:
  somatic spike   = upward crossing of 0 mV in V_TRACE (recomputed; count compared with SPIKES.tsv)
  dendritic spike = upward crossing of 0 mV in V_DEND_TRACE
  quiescent gap   = somatic ISI >= 1000 ms
  firing phase    = first spike after a gap to the start of the next gap
  tonic segment   = phase start to the first dendritic spike in the phase; burst segment = the remainder
  complete cycle  = phase with tonic >= 1 s, burst with >= 1 dendritic spike, followed by a gap >= 1 s
Sorter label: pinned ${REPO_ROOT}/05_model_forrest/p4_t2_regime.py (sha fdab1cb9...) with Bounds(v_quiescent_max=-60.0,
v_block_min=-30.0), applied to the full 90 s trace and, for information, to t >= 30 s.
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/stageA_rescore'))
SORTER = Path(__file__).resolve().parent / 'p4_t2_regime.py'
SORTER_SHA = 'b3714ca27132582b01fbfcf31ab9d4e25e95764af8292ec82027cab1cca71885'
BOUNDS = {'v_quiescent_max': -60.0, 'v_block_min': -30.0}
RUNS = {
    'run1_vinit_m70': {
        'dir': _release_paths.expand('${DATA_ROOT}/results/p4_t2_stageA_run_v02'),
        'V_TRACE.tsv': 'a84fe9ef8ced6a196cca56f14cafdec0b76d09226bf9644538c0de562cb205b8',
        'V_DEND_TRACE.tsv': '39af906511aa4dd9e4a9df9ba00b7bfc54de5d78fcf901a6e6fa838120bf26e1',
        'SPIKES.tsv': '2d7473a82c4f686ab1c20172f048fb2a373e10482d701a021e01c1e6b01548b5',
    },
    'run2_vinit_m65': {
        'dir': _release_paths.expand('${DATA_ROOT}/results/p4_t2_stageA_run_v03_vinit65'),
        'V_TRACE.tsv': '0890b56cc9475fa3a49f807c1d6638d2d5a4bb80c8318a15c414cef45e763d5e',
        'V_DEND_TRACE.tsv': 'e923801422e54e63cbe1132417508a8f12ccebb9d244adc80d9849928b9d9353',
        'SPIKES.tsv': 'b147fe22eced20b0ae9567cf98da08d71d6ef768145881094ae5246c1aa830b5',
    },
}
THRESHOLD_MV = 0.0
GAP_MS = 1000.0
TONIC_MIN_MS = 1000.0
SETTLE_INFO_MS = 30000.0


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f'not JSON serialisable: {type(o).__name__}')


def load_checked(path, expected):
    data = path.read_bytes()
    actual = digest(data)
    if actual != expected:
        raise SystemExit(f'hash mismatch for {path}: {actual} != {expected}')
    return data, actual


def read_trace(data):
    rows = np.loadtxt(data.decode().splitlines(), delimiter='\t', skiprows=1)
    header = data.split(b'\n', 1)[0].decode()
    if header != 't_ms\tv_mV':
        raise SystemExit(f'unexpected trace header {header!r}')
    return rows[:, 0], rows[:, 1]


def read_spikes(data):
    header, _, body = data.decode().partition('\n')
    if header != 'spike_ms':
        raise SystemExit(f'unexpected SPIKES header {header!r}')
    return np.array([float(x) for x in body.split()])


def upward_crossings(t, v, threshold=THRESHOLD_MV):
    """Times (linear interpolation) where v goes from below threshold to at or above it."""
    idx = np.flatnonzero((v[:-1] < threshold) & (v[1:] >= threshold))
    return t[idx] + (threshold - v[idx]) / (v[idx + 1] - v[idx]) * (t[idx + 1] - t[idx])


def mean_isi(times):
    return float(np.diff(times).mean()) if len(times) >= 2 else None


def score(t, v_soma, v_dend, spikes_file):
    soma = upward_crossings(t, v_soma)
    dend = upward_crossings(t, v_dend)
    isi = np.diff(soma)
    gap_idx = np.flatnonzero(isi >= GAP_MS)  # gap k lies between soma[gap_idx[k]] and soma[gap_idx[k]+1]
    gaps = [{'index': int(i), 'start_ms': float(soma[i]), 'end_ms': float(soma[i + 1]),
             'duration_ms': float(isi[i])} for i in gap_idx]
    periods = [float(gaps[k + 1]['start_ms'] - gaps[k]['start_ms']) for k in range(len(gaps) - 1)]
    # Phases: spikes strictly between consecutive gaps. Leading phase (trace start to first gap) and
    # trailing phase (last gap to trace end) are reported but flagged incomplete.
    bounds = [-1] + [int(i) for i in gap_idx] + [len(soma) - 1]
    phases = []
    for k in range(len(bounds) - 1):
        first, last = bounds[k] + 1, bounds[k + 1]
        if last < first:
            continue
        ph = soma[first:last + 1]
        start, end = float(ph[0]), float(ph[-1])
        d_in = dend[(dend >= start) & (dend <= end)]
        complete = k >= 1 and k < len(gaps)  # preceded and followed by a gap
        if len(d_in):
            split = float(d_in[0])
            tonic = ph[ph < split]
            burst = ph[ph >= split]
            tonic_ms, burst_ms = split - start, end - split
        else:
            tonic, burst = ph, ph[:0]
            tonic_ms, burst_ms = end - start, 0.0
        following_gap = float(gaps[k]['duration_ms']) if k < len(gaps) else None
        phases.append({
            'phase': k, 'preceded_by_gap': k >= 1, 'followed_by_gap': k < len(gaps),
            'start_ms': start, 'end_ms': end, 'duration_ms': end - start,
            'tonic_ms': float(tonic_ms), 'burst_ms': float(burst_ms),
            'n_soma_total': int(len(ph)), 'n_soma_tonic': int(len(tonic)), 'n_soma_burst': int(len(burst)),
            'mean_isi_tonic_ms': mean_isi(tonic), 'mean_isi_burst_ms': mean_isi(burst),
            'n_dend': int(len(d_in)), 'following_gap_ms': following_gap,
            'complete_cycle': bool(complete and tonic_ms >= TONIC_MIN_MS and len(d_in) >= 1
                                   and following_gap is not None and following_gap >= GAP_MS),
        })
    return {
        'trace_start_ms': float(t[0]), 'trace_end_ms': float(t[-1]), 'n_samples': int(len(t)),
        'dt_ms': float(np.median(np.diff(t))),
        'n_soma_spikes_recomputed': int(len(soma)), 'n_spikes_file': int(len(spikes_file)),
        'spike_count_difference': int(len(soma) - len(spikes_file)),
        'max_abs_time_diff_ms_if_same_count': (float(np.max(np.abs(soma - spikes_file)))
                                               if len(soma) == len(spikes_file) else None),
        'n_dend_spikes': int(len(dend)),
        'n_dend_outside_phases': int(len(dend) - sum(p['n_dend'] for p in phases)),
        'leading_silence_ms': float(soma[0] - t[0]) if len(soma) else None,
        'trailing_silence_ms': float(t[-1] - soma[-1]) if len(soma) else None,
        'n_gaps': len(gaps), 'gaps': gaps, 'gap_to_gap_periods_ms': periods,
        'phases': phases, 'n_complete_cycles': int(sum(p['complete_cycle'] for p in phases)),
    }


def sorter_labels(t, v):
    sorter_bytes = SORTER.read_bytes()
    if digest(sorter_bytes) != SORTER_SHA:
        raise SystemExit('sorter hash mismatch')
    sys.path.insert(0, str(SORTER.parent))
    import p4_t2_regime as sorter  # noqa: E402
    b = sorter.Bounds(**BOUNDS)
    out = {}
    for name, mask in (('full_window', t >= t[0]), (f'from_{int(SETTLE_INFO_MS)}ms', t >= SETTLE_INFO_MS)):
        out[name] = sorter.classify(list(t[mask]), list(v[mask]), b)  # epochs kept in RESCORE.json
    return {'sorter': str(SORTER.relative_to(ROOT)), 'sorter_sha256': SORTER_SHA, 'bounds': BOUNDS,
            'labels_available': list(sorter.LABELS), 'results': out}


def main():
    OUT.mkdir()  # exclusive: fails if the folder exists
    report = {'script': str(Path(__file__).relative_to(ROOT)), 'script_sha256': digest(Path(__file__).read_bytes()),
              'spec': 'pre-registered definitions, 2026-10-12', 'threshold_mV': THRESHOLD_MV,
              'gap_ms': GAP_MS, 'tonic_min_ms': TONIC_MIN_MS, 'inputs': {}, 'runs': {}}
    lines = []
    for run, spec in RUNS.items():
        d = ROOT / spec['dir']
        soma_b, h1 = load_checked(d / 'V_TRACE.tsv', spec['V_TRACE.tsv'])
        dend_b, h2 = load_checked(d / 'V_DEND_TRACE.tsv', spec['V_DEND_TRACE.tsv'])
        spk_b, h3 = load_checked(d / 'SPIKES.tsv', spec['SPIKES.tsv'])
        report['inputs'][run] = {'dir': spec['dir'], 'V_TRACE.tsv': h1, 'V_DEND_TRACE.tsv': h2, 'SPIKES.tsv': h3}
        t, vs = read_trace(soma_b)
        td, vd = read_trace(dend_b)
        if not np.array_equal(t, td):
            raise SystemExit(f'{run}: soma and dendrite time axes differ')
        s = score(t, vs, vd, read_spikes(spk_b))
        s['sorter'] = sorter_labels(t, vs)
        report['runs'][run] = s
        lines.append(f'{run}  ({spec["dir"]})')
        lines.append(f'  samples {s["n_samples"]}  dt {s["dt_ms"]:.3f} ms  span {s["trace_start_ms"]:.1f}-{s["trace_end_ms"]:.1f} ms')
        lines.append(f'  somatic spikes recomputed {s["n_soma_spikes_recomputed"]}  SPIKES.tsv {s["n_spikes_file"]}  '
                     f'difference {s["spike_count_difference"]}  max |dt| {s["max_abs_time_diff_ms_if_same_count"]}')
        lines.append(f'  dendritic spikes {s["n_dend_spikes"]}  gaps(ISI>={GAP_MS:.0f} ms) {s["n_gaps"]}  '
                     f'leading silence {s["leading_silence_ms"]:.1f} ms  trailing silence {s["trailing_silence_ms"]:.1f} ms')
        lines.append('  gap durations ms: ' + ', '.join(f'{g["duration_ms"]:.1f}' for g in s['gaps']))
        lines.append('  gap-to-gap periods ms: ' + ', '.join(f'{p:.1f}' for p in s['gap_to_gap_periods_ms']))
        for p in s['phases']:
            flag = 'complete' if p['complete_cycle'] else ('edge' if not (p['preceded_by_gap'] and p['followed_by_gap']) else 'incomplete')
            lines.append(f'  phase {p["phase"]:2d} [{flag:10s}] {p["start_ms"]:9.1f}-{p["end_ms"]:9.1f}  tonic {p["tonic_ms"]:8.1f} ms '
                         f'(n {p["n_soma_tonic"]:4d}, ISI {p["mean_isi_tonic_ms"]})  burst {p["burst_ms"]:8.1f} ms '
                         f'(n {p["n_soma_burst"]:4d}, ISI {p["mean_isi_burst_ms"]})  dend {p["n_dend"]:3d}  next gap {p["following_gap_ms"]}')
        lines.append(f'  complete cycles (tonic>=1 s, >=1 dendritic spike, gap>=1 s): {s["n_complete_cycles"]}')
        for w, r in s['sorter']['results'].items():
            lines.append(f'  sorter {w}: {r["label"]} ({r["reason"]}; n_spikes {r.get("n_spikes")}, rate {r.get("rate_hz")}, '
                         f'gaps {r.get("n_gaps")}, modes {r.get("modes")})')
        lines.append('')
    with (OUT / 'RESCORE.json').open('x') as f:
        f.write(json.dumps(report, indent=2, default=json_default) + '\n')
    with (OUT / 'REPORT.txt').open('x') as f:
        f.write('\n'.join(lines) + '\n')
    manifest = ''.join(f'{digest(p.read_bytes())}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file())
    with (OUT / 'SHA256SUMS.txt').open('x') as f:
        f.write(manifest)
    print('\n'.join(lines))
    print(manifest)


if __name__ == '__main__':
    main()

#%% debug
