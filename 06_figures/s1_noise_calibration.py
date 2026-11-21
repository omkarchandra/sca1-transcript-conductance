
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Supplementary Figure S1, WT noise calibration (
# The drawing script of the accepted Figure_S1_v2 (final_figures_2.1_phase1_to_phase4/supplementary, 14 Oct) was not found in the repository or under
# This redraw takes its layout from S1_v2 (log CV axis, twin rate axis, target CV line 0.08, sigma-used line 0.0065 nA,
# legend under the axes) and its data from CALIBRATION.tsv (6af9d5ce; scorer v08 bf8b41a3, masoli_pop_CAL_45s).
#   content check: the 5 sigma, 5 rates and 5 CVs as plotted, against the old S1_v2 SVG (marker positions mapped back with its own axis ticks).
#   font rule: built at 17.5 cm, never rescaled; ticks 8 pt, axis labels and title 9 pt, legend 8 pt; min font read back into RECORD.json.
#   outputs: PNG 300 dpi, PDF, SVG (text kept as text), RECORD.json, SHA256SUMS into a NEW real_v01/.
# Run: python -I -B 06_figures/s1_noise_calibration.py
import hashlib, json, math, re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/s1_noise_calibration'))
CAL = ROOT / _release_paths.expand('${REPO_ROOT}/results/calibration/CALIBRATION.tsv')
CAL_SHA = '6af9d5ce1f3b39ec3b3ce4d4b03261e0f0b64a1afc529f5c5e65fa47076d9093'
OLD = Path(_release_paths.expand('${WORK_ROOT}/Figure_S1_v2.svg'))
TARGET_CV, SIGMA_USED = 0.08, 0.0065
CV_COLOR, RATE_COLOR, TARGET_COLOR = '#0b7285', '#e67700', '#c92a2a'
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 10.0
TICK, LABEL, TITLE, LEGEND = 8, 9, 9, 8
NS = '{http://www.w3.org/2000/svg}'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    assert sha(CAL) == CAL_SHA, 'calibration table changed'
    t = pd.read_csv(CAL, sep='\t').sort_values('noise_sigma_nA').reset_index(drop=True)
    assert len(t) == 5 and t.condition.eq('WT').all() and t.run_class.eq('firing').all() and t.noise_sigma_nA.is_unique
    return t


def old_points():
    """Plotted (sigma, CV) and (sigma, rate) of the old S1_v2, from its marker positions and its own linear x / log-CV / linear rate ticks."""
    root = ET.parse(OLD).getroot(); axes = [g for g in root.iter(NS + 'g') if re.fullmatch(r'axes_\d+', g.get('id', ''))]

    def ticks(ax, kind):
        out = []
        for g in (g for g in ax.iter(NS + 'g') if re.fullmatch(kind + r'tick_\d+', g.get('id', ''))):
            u = list(g.iter(NS + 'use')); tx = ''.join(''.join(t.itertext()) for t in g.iter(NS + 'text')).strip()
            if u and tx: out.append((float(u[0].get('x' if kind == 'x' else 'y')), tx))
        return out

    def lin(tk):
        (p1, v1), (p2, v2) = (tk[0][0], float(tk[0][1].replace('−', '-'))), (tk[-1][0], float(tk[-1][1].replace('−', '-')))
        return lambda p: v1 + (p - p1) * (v2 - v1) / (p2 - p1)
    ax_cv, ax_rate = axes[0], axes[1]
    fx = lin(ticks(ax_cv, 'x'))
    # log CV axis: major tick labels are mathtext 10^k, rendered as separate text pieces; use their y with exponents -2, -1, 0 in order
    ymaj = [p for p, _ in ticks(ax_cv, 'y')]
    ys = sorted(set(round(p, 4) for p in ymaj))[::-1]  # pixel y grows downward: largest y = smallest value
    assert len(ys) >= 3
    e = {ys[0]: -2, ys[1]: -1, ys[2]: 0}
    (pa, ea), (pb, eb) = (ys[0], -2), (ys[2], 0)
    fcv = lambda p: 10 ** (ea + (p - pa) * (eb - ea) / (pb - pa))
    frate = lin(ticks(ax_rate, 'y'))
    def marks(ax):
        pts = []
        for line in (g for g in ax.iter(NS + 'g') if g.get('id', '').startswith('line2d_')):
            us = [u for u in line.iter(NS + 'use') if u.get('x') is not None]
            if len(us) == 5: pts.append([(float(u.get('x')), float(u.get('y'))) for u in us])
        return pts
    cv = marks(ax_cv); rate = marks(ax_rate)
    assert len(cv) == 1 and len(rate) == 1, (len(cv), len(rate))
    return [fx(x) for x, _ in cv[0]], [fcv(y) for _, y in cv[0]], [fx(x) for x, _ in rate[0]], [frate(y) for _, y in rate[0]]


def draw(t):
    fig, ax = plt.subplots(figsize=(WIDTH_CM * CM, HEIGHT_CM * CM), layout='constrained')
    ax.plot(t.noise_sigma_nA, t.isi_cv, 'o-', color=CV_COLOR, lw=1.4, ms=5)
    ax.axhline(TARGET_CV, ls='--', color=TARGET_COLOR, lw=1.2); ax.axvline(SIGMA_USED, ls=':', color='#555555', lw=1.2)
    ax.set_yscale('log'); ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}')); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())  # plain labels: mathtext exponents would print below 8 pt
    ax.set_xlabel('Noise σ (nA)', fontsize=LABEL); ax.set_ylabel('WT ISI CV (log scale)', fontsize=LABEL); ax.tick_params(labelsize=TICK)
    ax2 = ax.twinx(); ax2.plot(t.noise_sigma_nA, t.rate_hz, 's:', color=RATE_COLOR, lw=1.4, ms=5); ax2.set_ylabel('WT rate (Hz)', fontsize=LABEL); ax2.tick_params(labelsize=TICK)
    ax.set_title('WT noise calibration', fontsize=TITLE)
    ax.set_xticks(np.arange(0, 0.0401, 0.005)); ax.set_yticks([0.01, 0.1, 1]); ax2.set_yticks(np.arange(32, 44.1, 2))  # the ticks of S1_v2, all inside the view
    handles = [Line2D([], [], marker='o', ls='-', color=CV_COLOR, label='WT ISI CV (left axis)'), Line2D([], [], ls='--', color=TARGET_COLOR, label=f'Target CV ({TARGET_CV:g})'),
               Line2D([], [], ls=':', color='#555555', label=f'σ used ({SIGMA_USED:g} nA)'), Line2D([], [], marker='s', ls=':', color=RATE_COLOR, label='WT rate (right axis)')]
    fig.legend(handles=handles, loc='outside lower center', ncol=4, fontsize=LEGEND, frameon=False)
    return fig


def main():
    t = load()
    osx, ocv, orx, orate = old_points()
    check = dict(sigma_new=t.noise_sigma_nA.tolist(), sigma_old_cv_markers=osx, sigma_old_rate_markers=orx, rate_new=t.rate_hz.tolist(), rate_old=orate,
                 cv_new=t.isi_cv.tolist(), cv_old=ocv)
    rel = lambda a, b: float(np.max(np.abs(np.array(a) - np.array(b)) / np.maximum(1e-12, np.abs(np.array(b)))))
    check['max_abs_diff_sigma'] = float(max(np.max(np.abs(np.array(osx) - t.noise_sigma_nA)), np.max(np.abs(np.array(orx) - t.noise_sigma_nA))))
    check['max_abs_diff_rate_hz'] = float(np.max(np.abs(np.array(orate) - t.rate_hz))); check['max_rel_diff_cv'] = rel(ocv, t.isi_cv)
    check['n_compared'] = 15; check['n_equal'] = int((np.abs(np.array(osx) - t.noise_sigma_nA) < 1e-6).sum() + (np.abs(np.array(orate) - t.rate_hz) < 1e-3).sum() + (np.abs(np.array(ocv) - t.isi_cv) / t.isi_cv < 1e-4).sum())
    if check['n_equal'] != 15:
        raise SystemExit(f'STOP: content check differs: {json.dumps(check)}')
    fig = draw(t); fig.canvas.draw()
    texts = [x for x in fig.findobj(matplotlib.text.Text) if x.get_visible() and x.get_text().strip()]
    sizes = [x.get_fontsize() for x in texts]
    r = fig.canvas.get_renderer(); fw, fh = fig.bbox.width, fig.bbox.height
    overflow = [x.get_text() for x in texts if (e := x.get_window_extent(r)).x0 < -0.5 or e.x1 > fw + 0.5 or e.y0 < -0.5 or e.y1 > fh + 0.5]
    assert not overflow, ('text runs past the figure edge', overflow)
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_S1_noise_calibration.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_S1_noise_calibration.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_S1_noise_calibration.svg', metadata={'Date': None}); plt.close(fig)
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    rec = dict(task='Supplementary S1 rebuild', script_sha256=sha(Path(__file__)),
               note='The drawing script of Figure_S1_v2 was not found; this redraw takes its layout from S1_v2 and its data from CALIBRATION.tsv.',
               calibration_table=dict(path=str(CAL), sha256=CAL_SHA), old_figure=dict(path=str(OLD), sha256=sha(OLD)), target_cv=TARGET_CV, sigma_used_nA=SIGMA_USED,
               content_check=check, smallest_font_pt=min(sizes), font_sizes_used=sorted(set(sizes)), width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), dpi=300,
               matplotlib=matplotlib.__version__, outputs=outs, finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f: f.write(json.dumps(rec, indent=1) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS': f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(dict(content_check=check, smallest_font_pt=min(sizes), font_sizes=sorted(set(sizes)), width_cm=w_cm, height_cm=h_cm, **outs), indent=1))


if __name__ == '__main__':
    main()

#%% debug
