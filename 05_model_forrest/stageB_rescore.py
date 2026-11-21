""" blind re-score of the Track 2 Stage B runs (pre-specified).

python -I -B 05_model_forrest/stageB_rescore.py
Reads only V_TRACE.tsv, V_DEND_TRACE.tsv and SPIKES.tsv of the sixteen pinned Stage B folders (every input
hashed before parsing, hashes must equal the pinned values); never RECORD.json, SHA256SUMS, the summaries or
anything under docs/analysis_loop/. Definitions are written from review's spec, not from the code. Writes only
under ${DATA_ROOT}/analysis/stageB_rescore/ with exclusive create. No simulation, nothing deleted.

Definitions:
  window          = t >= 60,000 ms (first 60 s settle excluded from every number)
  somatic spike   = upward crossing of 0 mV in V_TRACE within the window (compared with SPIKES.tsv entries in window)
  dendritic spike = upward crossing of 0 mV in V_DEND_TRACE within the window (count only)
  gap             = interval between consecutive somatic spikes >= 1,000 ms; leading/trailing edges flagged
                    for every gap and edge: duration, mean/min/max V (soma), mean V (dendrite); no spike -> whole window;
                    stats of a gap or trailing edge start POST_SPIKE_MARGIN_MS after the opening spike (v02, review )
  period          = gap start to next gap start
  fraction above -40 mV = share of window V_TRACE samples with v > -40
  sorter          = pinned p4_t2_regime.py classify(list(t), list(v), Bounds(-60.0, -30.0)) on the window
  class (<AUTHOR> 00:52) = plateau if DEPOLARIZATION_BLOCK; silent if QUIESCENT; firing if >= 1 somatic spike;
                    else UNCLASSIFIED with the raw label
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np

ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/stageB_rescore'))
RUN_PREFIX = _release_paths.expand('${DATA_ROOT}/results/p4_t2_stageB_run_')
SORTER = Path(__file__).resolve().parent / 'p4_t2_regime.py'
SORTER_SHA = 'b3714ca27132582b01fbfcf31ab9d4e25e95764af8292ec82027cab1cca71885'
BOUNDS = {'v_quiescent_max': -60.0, 'v_block_min': -30.0}
WINDOW_START_MS = 60000.0
THRESHOLD_MV = 0.0
GAP_MS = 1000.0
ABOVE_MV = -40.0
POST_SPIKE_MARGIN_MS = 10.0  # interval stats after an opening spike skip its peak
FILES = ('V_TRACE.tsv', 'V_DEND_TRACE.tsv', 'SPIKES.tsv')
# Pinned input hashes : condition -> (V_TRACE, V_DEND_TRACE, SPIKES)
INPUTS = {
    'KI_18w_high': ('8faec7cfff57258a424a4aac4159cc51fab335a688539b4aff29ea3ca119fc92', '0a06e8b36a9cf6d4116045d94ab016a41300d42f35ccb71b2078c624c3f02512', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_18w_low': ('efd72538ed67c78cfd857b172296dfe6e1ca5105d7f12eb70ced3b30b42868bf', '16ce4985c26de83637ceabdbef07d12078feabeb70f86d09d0e7af21f1c14ec9', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_18w_mid': ('558fd286aced019b430afd29535d3777ae0ae5bab6b947ef0ea113a9cb2c1bfe', '4edb97c19d1cd094f8e4b599289923dfecdc7438b393f98d429f425d342650f7', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_24w_high': ('b5d26de08f191ee5bdd14f893b1df4c300f7f69d82032138eb306cb1280caf12', '4d3b1dd5aa4b04c6354d2be126f50a311b4f95e8fa4e3982bea033c492db7325', '701315996e7c9b0124b04fe0aabe0cfd058fb5ed95b6f58f27b2deba0266b68b'),
    'KI_24w_low': ('535ed4b5e815f17bde63975eda2e538983b84559aaddeef2c0256f493c25264f', '522f2d4856e828345c6ba0bc5cd7869c06a7775b87166ee52c74877ee44d69db', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_24w_mid': ('bf6cb469149aed8382a675e668f66736f5e4003c05b14ee56c6d9d4d40ff802a', '7fb8b14e7a3d922d268680a07ecb17721bf5e96dd16321b417175b905b57d03d', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_30w_high': ('185eb77e36192aaf36f55c0ee8a7c3b4136570027ca5894229185f5ac3499467', '48faef52f369a49a38061ef500a9291ed9348da3cc6cae730888db76cc5f5eae', '42edf84e754b8350a2e5199928d2b9325b6b00c27f917434ca0c244f8ed3e867'),
    'KI_30w_low': ('c48188fcb56c0fe15ca5e75615605e08c285a51116b8562d64183f0d160da876', '2b9947a41ba63acf166c60c09696c0051d089f68a131e3308fd912b32466fd24', 'ac99ce282a212eef68bbf50a60269aecc8d3dcc9fb51c6c8bb4a7c7c0dff8e26'),
    'KI_30w_mid': ('b55182160edbc8fb22aafc10058f9cd5c4838d912860fe2e302e9620a995f636', 'de413d4879da55efe62d27ac61bfe06a8cfbaadbb0dd00d37dea844af5811fb5', 'cb9f8315dd4bd96e4da6eb3c1baa03c3485d7a196cbfb85c4dff81afe358edb8'),
    'KI_5w_high': ('76ed5c561c6c585f96284e731a1babe43b1a9fc3afd7091143c8fd97217ff519', 'b83b1ed00bd217b113f7dc284988fe7dd3502bd2ac46625e8315dee4467b698e', '1f4c54f3cb4df773712e0db1c67241d0627d9e1bfce61a389828707f56fd9e9e'),
    'KI_5w_low': ('8094616f6db8424cfce1679792df8afaaba1dc981cf532d32d395a2ca3386291', '6cf7c005288c4e0f4140bf370d306f0f2ad1717c93646b5f01cccc81c2a0c1a3', 'c432277f5fc0ed4ee14b16ed8d4f00a6671838009e0b79c75bc26594597e750f'),
    'KI_5w_mid': ('3ee13e52810deb0011dbcbcfc640f78605a0e04ce81eef6a17db353d8d849b02', 'b3a29005458d0f3e1dcdb33667060f74bb1d380c7d1edd18244edc0df866661b', '6b6738f5aa6aecf8669c3500c44a77976905f9e8dee0a74bc294edbb6c677e4a'),
    'WT_18w_reference': ('a6f7052ab2cca8ce3f1c15031954e0ba807852fe7f070594c5e9f3833f3d42da', '166e60131285f30c4d8cb87fccad430d670ac2bff654b33e7a965c1a4b2d9780', '4692f593e2b83350dfc5b35ccbcb55f7fdde2ed4cc77b025da28d5005520ad9d'),
    'WT_24w_reference': ('a6f7052ab2cca8ce3f1c15031954e0ba807852fe7f070594c5e9f3833f3d42da', '166e60131285f30c4d8cb87fccad430d670ac2bff654b33e7a965c1a4b2d9780', '4692f593e2b83350dfc5b35ccbcb55f7fdde2ed4cc77b025da28d5005520ad9d'),
    'WT_30w_reference': ('a6f7052ab2cca8ce3f1c15031954e0ba807852fe7f070594c5e9f3833f3d42da', '166e60131285f30c4d8cb87fccad430d670ac2bff654b33e7a965c1a4b2d9780', '4692f593e2b83350dfc5b35ccbcb55f7fdde2ed4cc77b025da28d5005520ad9d'),
    'WT_5w_reference': ('a6f7052ab2cca8ce3f1c15031954e0ba807852fe7f070594c5e9f3833f3d42da', '166e60131285f30c4d8cb87fccad430d670ac2bff654b33e7a965c1a4b2d9780', '4692f593e2b83350dfc5b35ccbcb55f7fdde2ed4cc77b025da28d5005520ad9d'),
}
ORDER = ['WT_5w_reference', 'WT_18w_reference', 'WT_24w_reference', 'WT_30w_reference'] + [
    f'KI_{age}_{level}' for age in ('5w', '18w', '24w', '30w') for level in ('low', 'mid', 'high')]


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
    return data


def read_trace(data):
    header, _, body = data.partition(b'\n')
    if header != b't_ms\tv_mV':
        raise SystemExit(f'unexpected trace header {header!r}')
    rows = np.loadtxt(io.BytesIO(body), delimiter='\t', ndmin=2)
    return rows[:, 0], rows[:, 1]


def read_spikes(data):
    header, _, body = data.decode().partition('\n')
    if header != 'spike_ms':
        raise SystemExit(f'unexpected SPIKES header {header!r}')
    return np.array([float(x) for x in body.split()], dtype=float)


def upward_crossings(t, v, threshold=THRESHOLD_MV):
    """Times (linear interpolation) where v goes from below threshold to at or above it."""
    idx = np.flatnonzero((v[:-1] < threshold) & (v[1:] >= threshold))
    return t[idx] + (threshold - v[idx]) / (v[idx + 1] - v[idx]) * (t[idx + 1] - t[idx])


def interval_stats(t, vs, vd, start, end, kind, after_spike):
    stats_from = start + (POST_SPIKE_MARGIN_MS if after_spike else 0.0)
    sel = (t >= stats_from) & (t <= end)
    return {'kind': kind, 'start_ms': float(start), 'end_ms': float(end), 'duration_ms': float(end - start),
            'stats_from_ms': float(stats_from),
            'n_samples': int(sel.sum()),
            'mean_v_soma': float(vs[sel].mean()) if sel.any() else None,
            'min_v_soma': float(vs[sel].min()) if sel.any() else None,
            'max_v_soma': float(vs[sel].max()) if sel.any() else None,
            'mean_v_dend': float(vd[sel].mean()) if sel.any() else None}


def score(t, vs, vd, spikes_file):
    full_n, dt = len(t), float(np.median(np.diff(t)))
    w = t >= WINDOW_START_MS
    t, vs, vd = t[w], vs[w], vd[w]
    t0, t1 = float(t[0]), float(t[-1])
    soma = upward_crossings(t, vs)
    dend = upward_crossings(t, vd)
    file_in_window = spikes_file[(spikes_file >= t0) & (spikes_file <= t1)]
    intervals = []
    if len(soma) == 0:
        intervals.append(interval_stats(t, vs, vd, t0, t1, 'whole_window_no_spike', False))
    else:
        intervals.append(interval_stats(t, vs, vd, t0, float(soma[0]), 'leading_edge', False))
        isi = np.diff(soma)
        for i in np.flatnonzero(isi >= GAP_MS):
            intervals.append(interval_stats(t, vs, vd, float(soma[i]), float(soma[i + 1]), 'gap', True))
        intervals.append(interval_stats(t, vs, vd, float(soma[-1]), t1, 'trailing_edge', True))
    gaps = [x for x in intervals if x['kind'] == 'gap']
    periods = [float(gaps[k + 1]['start_ms'] - gaps[k]['start_ms']) for k in range(len(gaps) - 1)]
    return {
        'n_samples_full_trace': int(full_n), 'dt_ms': dt,
        'window_first_sample_ms': t0, 'window_last_sample_ms': t1, 'n_samples_window': int(len(t)),
        'n_soma_spikes_recomputed': int(len(soma)), 'n_spikes_file_in_window': int(len(file_in_window)),
        'n_spikes_file_total': int(len(spikes_file)),
        'spike_count_difference': int(len(soma) - len(file_in_window)),
        'max_abs_time_diff_ms_if_same_count': (float(np.max(np.abs(soma - file_in_window)))
                                               if len(soma) and len(soma) == len(file_in_window) else None),
        'n_dend_spikes': int(len(dend)),
        'mean_v': float(vs.mean()), 'min_v': float(vs.min()), 'max_v': float(vs.max()),
        'fraction_above_m40': float((vs > ABOVE_MV).mean()),
        'n_gaps': len(gaps), 'longest_gap_ms': max((g['duration_ms'] for g in gaps), default=None),
        'intervals': intervals, 'periods_ms': periods,
    }, t, vs


def load_sorter():
    if digest(SORTER.read_bytes()) != SORTER_SHA:
        raise SystemExit('sorter hash mismatch')
    sys.path.insert(0, str(SORTER.parent))
    import p4_t2_regime as sorter  # noqa: E402
    return sorter


def stage_b_class(label, n_spikes):
    if label == 'DEPOLARIZATION_BLOCK':
        return 'plateau'
    if label == 'QUIESCENT':
        return 'silent'
    if n_spikes >= 1:
        return 'firing'
    return f'UNCLASSIFIED ({label})'


def fmt(x, nd=2):
    return 'None' if x is None else f'{x:.{nd}f}'


def main():
    sorter = load_sorter()  # all inputs are hashed, loaded and scored before the output folder is created
    bounds = sorter.Bounds(**BOUNDS)
    report = {'script': str(Path(__file__).relative_to(ROOT)), 'script_sha256': digest(Path(__file__).read_bytes()),
              'spec': 'pre-registered definitions, 2026-10-12', 'window_start_ms': WINDOW_START_MS,
              'threshold_mV': THRESHOLD_MV, 'gap_ms': GAP_MS, 'above_mV': ABOVE_MV,
              'post_spike_margin_ms': POST_SPIKE_MARGIN_MS,
              'sorter': str(SORTER.relative_to(ROOT)), 'sorter_sha256_computed': digest(SORTER.read_bytes()),
              'sorter_bounds': BOUNDS, 'sorter_labels_available': list(sorter.LABELS),
              'class_rule': 'plateau=DEPOLARIZATION_BLOCK; silent=QUIESCENT; firing=>=1 somatic spike; else UNCLASSIFIED',
              'inputs': {}, 'conditions': {}}
    lines, rows, manifest = [], [], []
    for cond in ORDER:
        folder = ROOT / (RUN_PREFIX + cond)
        datas = {}
        for name, expected in zip(FILES, INPUTS[cond]):
            datas[name] = load_checked(folder / name, expected)
            manifest.append(f'{expected}  {RUN_PREFIX}{cond}/{name}')
        report['inputs'][cond] = dict(zip(FILES, INPUTS[cond]))
        t, vs = read_trace(datas['V_TRACE.tsv'])
        td, vd = read_trace(datas['V_DEND_TRACE.tsv'])
        if not np.array_equal(t, td):
            raise SystemExit(f'{cond}: soma and dendrite time axes differ')
        s, tw, vw = score(t, vs, vd, read_spikes(datas['SPIKES.tsv']))
        res = sorter.classify(list(tw), list(vw), bounds)
        s['sorter'] = res
        s['stage_b_class'] = stage_b_class(res['label'], s['n_soma_spikes_recomputed'])
        report['conditions'][cond] = s
        lines.append(f'{cond}  ({RUN_PREFIX}{cond})')
        lines.append(f'  full trace {s["n_samples_full_trace"]} samples, dt {s["dt_ms"]:.4f} ms; window {s["window_first_sample_ms"]:.1f}-'
                     f'{s["window_last_sample_ms"]:.1f} ms ({s["n_samples_window"]} samples)')
        lines.append(f'  somatic spikes recomputed {s["n_soma_spikes_recomputed"]}  SPIKES.tsv in window {s["n_spikes_file_in_window"]} '
                     f'(file total {s["n_spikes_file_total"]})  difference {s["spike_count_difference"]}  max |dt| {s["max_abs_time_diff_ms_if_same_count"]}')
        lines.append(f'  dendritic spikes {s["n_dend_spikes"]}  mean V {s["mean_v"]:.2f}  min V {s["min_v"]:.2f}  max V {s["max_v"]:.2f}  '
                     f'fraction above -40 mV {s["fraction_above_m40"]:.4f}')
        lines.append(f'  gaps (ISI >= {GAP_MS:.0f} ms) {s["n_gaps"]}  longest gap {fmt(s["longest_gap_ms"], 1)} ms  periods ms: '
                     + ', '.join(f'{p:.1f}' for p in s['periods_ms']))
        for x in s['intervals']:
            lines.append(f'    {x["kind"]:22s} {x["start_ms"]:9.1f}-{x["end_ms"]:9.1f}  dur {x["duration_ms"]:8.1f} ms  stats from {x["stats_from_ms"]:9.1f}  soma mean {fmt(x["mean_v_soma"])} '
                         f'min {fmt(x["min_v_soma"])} max {fmt(x["max_v_soma"])}  dend mean {fmt(x["mean_v_dend"])}')
        lines.append(f'  sorter: {res["label"]} ({res["reason"]})  n_spikes {res.get("n_spikes")}  rate {res.get("rate_hz")}  '
                     f'gaps {res.get("n_gaps")}  modes {res.get("modes")}')
        lines.append(f'  Stage B class: {s["stage_b_class"]}')
        lines.append('')
        rows.append((cond, s['n_soma_spikes_recomputed'], s['n_dend_spikes'], s['n_gaps'], fmt(s['longest_gap_ms'], 1),
                     f'{s["fraction_above_m40"]:.4f}', f'{s["mean_v"]:.2f}', res['label'], s['stage_b_class']))
    header = ('condition', 'spikes', 'dendritic', 'gaps>=1s', 'longest_gap_ms', 'frac_above_-40', 'mean_V', 'sorter_label', 'class')
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    lines.append('TABLE (window t >= 60 s)')
    for r in [header] + rows:
        lines.append('  ' + '  '.join(str(c).ljust(w) for c, w in zip(r, widths)))
    OUT.mkdir()  # exclusive: fails if the folder exists
    with (OUT / 'RESCORE.json').open('x') as f:
        f.write(json.dumps(report, indent=2, default=json_default) + '\n')
    with (OUT / 'REPORT.txt').open('x') as f:
        f.write('\n'.join(lines) + '\n')
    for name in ('RESCORE.json', 'REPORT.txt'):
        manifest.append(f'{digest((OUT / name).read_bytes())}  {name}')
    with (OUT / 'SHA256SUMS').open('x') as f:
        f.write('\n'.join(manifest) + '\n')
    print('\n'.join(lines))
    print('\n'.join(manifest[-2:]))


if __name__ == '__main__':
    main()

#%% debug
