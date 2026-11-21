
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Figure 3 (Masoli model), rebuilt to the font rule (
# Numbers: the accepted E4 panel scripts' own read_numbers() (fig5a_values.py d1cebe2a, fig5b_values.py 998a83b3; both pin the
# 840-run scoring SUMMARY abf0e8cf, PER_MEMBER 148890d0, PER_RUN b8b9f2a8 and recheck medians and counts against SUMMARY). Titles, axis labels,
# x labels, legend, colours, positions, jitter, threshold ticks, medians, the 0.75-1.25 band and the n/14 counts are those of the accepted panels.
#   layout: panel A above panel B, each at the full 17.5 cm width (the side-by-side Figure_3_v3 put each panel at 8.4 cm).
#   font rule: built at 17.5 cm, never rescaled; ticks and counts 8 pt, axis labels and titles 9 pt, legend 8 pt, panel letters 12 pt bold;
#           min font read back into RECORD.json.  outputs: PNG 300 dpi, PDF, SVG (text kept as text), RECORD.json, SHA256SUMS, NEW real_v01/.
# Run: python -I -B 06_figures/fig5_masoli.py
import hashlib, json, runpy
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/fig5_masoli'))
E4 = ROOT / _release_paths.expand('${REPO_ROOT}/06_figures')
SRC = {'A': (E4 / 'fig5a_values.py', '3d8e486990ad24c8a5f93b2a6730d70406516fe1bd8a8a8a6f46e2cb9dea339a'),
       'B': (E4 / 'fig5b_values.py', '9163ee6d6c990b523db18c8d1cbb27cd8aea21cfd93e0596aac3ae910ddb6087')}
SCORED = Path(_release_paths.expand('${REPO_ROOT}/results/e4'))
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 20.0
TICK, LABEL, TITLE, LEGEND, LETTER = 8, 9, 9, 8, 12


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    mods = {}
    for k, (p, want) in SRC.items():
        assert sha(p) == want, ('E4 source script changed', p)
        mods[k] = runpy.run_path(str(p), run_name='e4_source')
    a, ha = mods['A']['read_numbers'](SCORED); b, hb = mods['B']['read_numbers'](SCORED); assert ha == hb
    return mods, a, b, ha


def panel_a(ax, m, plotted):
    palette = {5: '#174e8c', 18: '#087f8c', 24: '#7b4ea3', 30: '#c35b20'}
    for age, lo, hi in ((18, 1.04, 3.96), (30, 8.04, 10.96)):
        ax.axvspan(lo, hi, color='#edf2f5', zorder=0)
    ax.axhline(0, color='#343a40', linewidth=.7, zorder=1)
    for entry in plotted:
        x = entry['x']
        for rank, point in enumerate(entry['points']):
            px = x + (rank - 6.5) * .025
            ax.plot([px - .012, px + .012], [point['cv_margin']] * 2, color='#777e85', linewidth=.7, zorder=2)
            ax.scatter(px, point['cv_increase'], s=14, color=palette[entry['age']], edgecolors='#22303a', linewidths=.3, zorder=3)
        assert entry['median_cv_increase'] == median(p['cv_increase'] for p in entry['points'])
        ax.plot([x - .24, x + .24], [entry['median_cv_increase']] * 2, color='#1c2530', linewidth=1.5, zorder=4)
        ax.text(x, 1.01, f"{entry['n_increased']}/{entry['n']}", ha='center', va='bottom', fontsize=TICK, transform=ax.get_xaxis_transform())
    ax.set_xticks(m['POSITIONS'], ['mid'] + ['low', 'mid', 'high'] * 3, fontsize=TICK)
    for age, x in ((5, 0), (18, 2.5), (24, 6), (30, 9.5)):
        ax.annotate(f'{age} W', (x, 0), xycoords=ax.get_xaxis_transform(), xytext=(0, -20), textcoords='offset points', ha='center', va='top', fontsize=TICK)
    ax.set_xlim(-.55, 11.05); ax.set_ylabel(m['Y_LABEL'], fontsize=LABEL); ax.tick_params(axis='y', labelsize=TICK)
    ax.set_title(m['TITLE'], fontsize=TITLE, pad=16)


def panel_b(ax, m, plotted):
    dark, light = '#163b66', '#f4a261'
    ax.axhspan(.75, 1.25, color='#ced4da', alpha=.35)
    ax.axhline(1, color='#495057', linewidth=.7)
    labels = []
    for i, entry in enumerate(plotted):
        color = light if entry['includes_Cacna1a'] else dark
        for rank, point in enumerate(entry['points']):
            ax.scatter(i + (rank - 6.5) * .012, point['rate_ratio'], s=14, facecolor=color, edgecolor='#1f2933', linewidth=.3, zorder=3)
        assert entry['median_rate_ratio'] == median(p['rate_ratio'] for p in entry['points'])
        ax.plot([i - .2, i + .2], [entry['median_rate_ratio']] * 2, color='#212529', linewidth=1.5, zorder=4)
        labels.append(m['LABELS'][entry['condition']])
    ax.set_xticks(range(len(plotted)), labels, fontsize=TICK, rotation=40, ha='right', rotation_mode='anchor')
    ax.set_xlim(-.55, len(plotted) - .45); ax.set_ylabel(m['Y_LABEL'], fontsize=LABEL); ax.tick_params(axis='y', labelsize=TICK)
    ax.set_title('KI firing-rate change and restorations (18 W)', fontsize=TITLE)
    return [Patch(color=light, label=m['LEGEND']), Patch(color='#ced4da', alpha=.35, label='reversal-rule rate tolerance (0.75–1.25)')]


def draw(mods, a, b):
    fig, (axa, axb) = plt.subplots(2, 1, figsize=(WIDTH_CM * CM, HEIGHT_CM * CM), height_ratios=[1, 1.15], layout='constrained')
    panel_a(axa, mods['A'], a); handles = panel_b(axb, mods['B'], b)
    fig.legend(handles=handles, loc='outside lower center', ncol=2, frameon=False, fontsize=LEGEND)
    for ax, letter in ((axa, 'A'), (axb, 'B')):
        ax.text(-0.09, 1.04, letter, transform=ax.transAxes, fontsize=LETTER, fontweight='bold', va='bottom', ha='left')
    return fig


def main():
    mods, a, b, hashes = load()
    fig = draw(mods, a, b); fig.canvas.draw()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = [t.get_fontsize() for t in texts]
    r = fig.canvas.get_renderer(); fw, fh = fig.bbox.width, fig.bbox.height
    overflow = [t.get_text() for t in texts if (e := t.get_window_extent(r)).x0 < -0.5 or e.x1 > fw + 0.5 or e.y0 < -0.5 or e.y1 > fh + 0.5]
    assert not overflow, ('text runs past the figure edge', overflow)
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_3_masoli.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_3_masoli.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_3_masoli.svg', metadata={'Date': None}); plt.close(fig)
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    rec = dict(task='Figure 3 (Masoli) font-rule rebuild', script_sha256=sha(Path(__file__)),
               e4_source_scripts={k: dict(path=str(p), sha256=h) for k, (p, h) in SRC.items()}, scored_dir=str(SCORED), input_sha256=hashes,
               medians_A={e['condition']: e['median_cv_increase'] for e in a}, counts_A={e['condition']: f"{e['n_increased']}/{e['n']}" for e in a},
               medians_B={e['condition']: e['median_rate_ratio'] for e in b},
               smallest_font_pt=min(sizes), font_sizes_used=sorted(set(sizes)), width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), dpi=300, matplotlib=matplotlib.__version__,
               outputs=outs, finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f: f.write(json.dumps(rec, indent=1) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS': f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(dict(smallest_font_pt=min(sizes), font_sizes=sorted(set(sizes)), width_cm=w_cm, height_cm=h_cm, **outs), indent=1))


if __name__ == '__main__':
    main()

#%% debug
