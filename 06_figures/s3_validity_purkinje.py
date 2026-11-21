
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Supplementary Figure S3, Purkinje-only (
# (the Phase 4 figure script validity_panel) with the 3 Purkinje scopes only (P, G(a), G(b)); the 3 Bergmann scopes dropped.
# Data: the verified Phase 4 closure bundle (job 9335707, SHARED_V02_COMPLETE_9335707.json 8a2feda9), Purkinje VALIDITY tables, same checks.
#   also writes LEGEND_COUNTS.json: path-position readouts passing rho >= 0.6 (of 3), DPT readouts passing (of 3), and per space the number of
#   the three WT 12-week mice between the WT 5- and 18-week medians (stored column WT_12w_mice_between_5_and_18).
#   font rule: built at 17.5 cm, never rescaled; ticks 8 pt, labels and titles 9 pt, legend and note 8 pt; min font read back into RECORD.json.
#   outputs: PNG 300 dpi, PDF, SVG (text kept as text), LEGEND_COUNTS.json, RECORD.json, SHA256SUMS into a NEW real_v01/.
# Run: python -I -B 06_figures/s3_validity_purkinje.py
import hashlib, io, json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/s3_validity_purkinje'))
BUNDLE = ROOT / _release_paths.expand('${DATA_ROOT}/phase4_shared_readout_rounding_review_v04/SHARED_V02_COMPLETE_9335707.json')
BUNDLE_SHA = '8a2feda9740166c4d00c96860c4e9f9db5f8f44dec475600791fad7c7aee7f78'
MANIFEST_SHA = 'e705ecaf14b1f297985dca25f2d6823a3d6e5da10575180c187359f0a0876e78'
SPACES = ('P', 'G_a', 'G_b'); POP = 'Purkinje cell'
POSITION_COLOR, DPT_COLOR = '#007F73', '#7851A9'
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 8.5
TICK, LABEL, TITLE, LEGEND, VALUE, VERDICT = 8, 9, 9, 8, 8, 9


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    data = BUNDLE.read_bytes(); assert hashlib.sha256(data).hexdigest() == BUNDLE_SHA, 'completion bundle changed'
    bundle = json.loads(data); manifest = bundle['SHA256SUMS.txt']['text'].encode(); assert hashlib.sha256(manifest).hexdigest() == MANIFEST_SHA
    pins = {name: h for h, name in (line.split(None, 1) for line in manifest.decode().splitlines())}
    val, used = {}, {}
    for space in SPACES:
        name = space + '_' + POP.replace(' ', '_') + '_VALIDITY.tsv'; entry = bundle['tables'][name]; content = entry['text'].encode()
        assert hashlib.sha256(content).hexdigest() == entry['sha256'] == pins[name], name
        v = pd.read_csv(io.BytesIO(content), sep='\t')
        assert v.space.eq(space).all() and v.population.eq(POP).all() and len(v) == 2 and set(v.readout) == {'position_t', 'DPT'} and v.threshold.eq(0.6).all()
        for r in v.itertuples(): assert r.informative == (r.spearman_rho_age_vs_crossfit_WT_mouse_median >= 0.6)
        val[space] = v.set_index('readout'); used[name] = pins[name]
    counts = dict(position_t_pass_of_3=int(sum(bool(val[s].loc['position_t', 'informative']) for s in SPACES)),
                  DPT_pass_of_3=int(sum(bool(val[s].loc['DPT', 'informative']) for s in SPACES)),
                  WT_12w_mice_between_5_and_18_medians_of_3={s: int(val[s].loc['position_t', 'WT_12w_mice_between_5_and_18']) for s in SPACES},
                  rho={s: {k: float(val[s].loc[k, 'spearman_rho_age_vs_crossfit_WT_mouse_median']) for k in ('position_t', 'DPT')} for s in SPACES})
    return val, used, counts


def draw(val):
    fig, (ax, sanity) = plt.subplots(1, 2, figsize=(WIDTH_CM * CM, HEIGHT_CM * CM), width_ratios=[3.0, 1.0], layout='constrained')
    ax.axvline(0.6, ls='--', color='#68717D', lw=1.0); ax.axvline(0, color='#E1E5EA', lw=0.8)
    for i, space in enumerate(SPACES):
        t = val[space]
        for name, color, marker, offset in (('position_t', POSITION_COLOR, 'o', -0.15), ('DPT', DPT_COLOR, 's', 0.15)):
            row = t.loc[name]; rho = float(row.spearman_rho_age_vs_crossfit_WT_mouse_median)
            ax.scatter(rho, i + offset, marker=marker, s=40, facecolor=color if row.informative else 'white', edgecolor=color, lw=1.3, zorder=3)
            right = rho < 0.78 and not 0.3 < rho < 0.6  # keep labels clear of the 0.6 threshold line
            ax.text(rho + (0.035 if right else -0.035), i + offset, f'{rho:.3f}', ha='left' if right else 'right', va='center', fontsize=VALUE, color=color)
        count = int(t.loc['position_t', 'WT_12w_mice_between_5_and_18'])
        sanity.text(0.5, i, f'{count} / 3 · ' + ('PASS' if count >= 2 else 'FAIL'), ha='center', va='center', fontsize=VERDICT, fontweight='bold', color='#9B3D43' if count < 2 else '#117733')
    ax.set_yticks(range(3), [s.replace('_a', '(a)').replace('_b', '(b)') + ' · Purkinje cell' for s in SPACES]); ax.tick_params(axis='y', length=0, pad=6, labelsize=LABEL)
    ax.set_ylim(2.6, -0.6); ax.set_xlim(-1.04, 1.04); ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.tick_params(axis='x', labelsize=TICK)
    ax.set_xlabel('Spearman ρ: age vs cross-fit WT mouse median', fontsize=LABEL)
    ax.set_title('Age ordering · required ρ ≥ 0.6', loc='left', fontsize=TITLE)
    ax.grid(axis='x', color='#EEF0F3', lw=0.6); ax.spines[['top', 'right', 'left']].set_visible(False)
    sanity.set_xlim(0, 1); sanity.set_ylim(ax.get_ylim()); sanity.set_title('WT 12 W position\nbetween 5 and 18 W\n(§27.3 needs ≥ 2 of 3)', fontsize=TITLE); sanity.set_axis_off()
    handles = [Line2D([], [], marker='o', ls='', color=POSITION_COLOR, label='path position'), Line2D([], [], marker='s', ls='', color=DPT_COLOR, label='DPT'),
               Line2D([], [], marker='o', ls='', color='#56616E', label='filled: ρ passes'), Line2D([], [], marker='o', ls='', mfc='white', mec='#56616E', label='hollow: ρ fails')]
    fig.legend(handles=handles, loc='outside lower center', ncol=4, fontsize=LEGEND, frameon=False,
               title='Right: number of the 3 WT 12-week mouse medians between the WT 5- and 18-week medians (stored check).', title_fontsize=LEGEND)
    return fig


def main():
    val, used, counts = load()
    fig = draw(val); fig.canvas.draw()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = [t.get_fontsize() for t in texts]
    r = fig.canvas.get_renderer(); fw = fig.bbox.width
    overflow = [t.get_text() for t in texts if t.get_window_extent(r).x0 < -0.5 or t.get_window_extent(r).x1 > fw + 0.5]
    assert not overflow, ('text runs past the figure edge', overflow)
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_S3_purkinje_only.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_S3_purkinje_only.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_S3_purkinje_only.svg', metadata={'Date': None}); plt.close(fig)
    with open(OUT / 'LEGEND_COUNTS.json', 'x') as f: f.write(json.dumps(counts, indent=1) + '\n')
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    rec = dict(task='Supplementary S3 Purkinje-only', script_sha256=sha(Path(__file__)), bundle=dict(path=str(BUNDLE), sha256=BUNDLE_SHA, manifest_sha256=MANIFEST_SHA),
               input_tables=used, legend_counts=counts, smallest_font_pt=min(sizes), font_sizes_used=sorted(set(sizes)), width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), dpi=300,
               matplotlib=matplotlib.__version__, outputs=outs, finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f: f.write(json.dumps(rec, indent=1) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS': f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(dict(counts=counts, smallest_font_pt=min(sizes), font_sizes=sorted(set(sizes)), width_cm=w_cm, height_cm=h_cm, **outs), indent=1))


if __name__ == '__main__':
    main()

#%% debug
