
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Supplementary Figure S4, Purkinje-only (
# S4 (final_figures_2.1_phase1_to_phase4/supplementary/Figure_S4.svg, 14 Oct) with only the 4 Purkinje rows, from HELDOUT_BY_PC.tsv (f32ed294).
#   columns as in the accepted S4: discovery KI - WT sign at 5/18/24/30 W, predicted sign, 12 W KI - WT (3 significant digits), sign match, complete separation.
#   gate:   sign matches 4 of 4 and complete separation 3 of 4 (PC5, PC6, PC30), else stop.
#   font rule: built at 17.5 cm, never rescaled; all text >= 9 pt except the 8 pt note line; title 10 pt; min font read back into RECORD.json.
#   outputs: PNG 300 dpi, PDF, SVG (text kept as text), RECORD.json, SHA256SUMS into a NEW real_v01/.
# Run: python -I -B 06_figures/s4_heldout_purkinje.py
import hashlib, json
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/s4_heldout_purkinje'))
HO = _release_paths.expand('${REPO_ROOT}/results/heldout/HELDOUT_BY_PC.tsv')
PINS = {HO: 'f32ed29456c1dbafae14369dcba9e358a90ae0d240b3d94a819f4e377ed73d5b'}
WT_COLOR, KI_COLOR, YES, NO = '#0072B2', '#D55E00', '#117733', '#888888'
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 7.3
BODY, HEAD, TITLE, NOTE, MARK = 9, 9, 10, 8, 11
MINUS = '−'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    for rel, h in PINS.items():
        assert sha(ROOT / rel) == h, ('pinned file changed', rel)
    ho = pd.read_csv(ROOT / HO, sep='\t'); ho = ho[ho.population == 'Purkinje_cell'].reset_index(drop=True); assert list(ho.PC) == [5, 6, 11, 30]
    n_match = int((ho.sign_match == 'yes').sum()); sep = [int(k) for k in ho.PC[ho.complete_separation_in_predicted_direction == 'yes']]
    if n_match != 4 or sep != [5, 6, 30]:
        raise SystemExit(f'STOP: held-out gate failed: sign matches {n_match}/4, complete separation {sep}')
    return ho, dict(sign_matches=n_match, complete_separation=sep)


def draw(ho):
    fig = plt.figure(figsize=(WIDTH_CM * CM, HEIGHT_CM * CM)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, WIDTH_CM); ax.set_ylim(0, HEIGHT_CM); ax.axis('off')
    # x positions in cm: (name, left, width)
    cols = [('PC', 0.3, 3.3), ('5 W', 3.7, 1.1), ('18 W', 4.9, 1.1), ('24 W', 6.1, 1.1), ('30 W', 7.3, 1.1), ('predicted', 8.7, 2.0), ('12 W KI ' + MINUS + ' WT', 11.0, 2.8), ('sign\nmatch', 14.0, 1.4), ('complete\nseparation', 15.5, 1.9)]
    top = HEIGHT_CM - 0.5
    ax.text(WIDTH_CM / 2, top, 'Held-out 12-week check (pre-declared rule), Purkinje cells', ha='center', va='center', fontsize=TITLE, fontweight='bold')
    ax.text(6.05, top - 0.8, 'discovery KI ' + MINUS + ' WT sign by age', ha='center', va='center', fontsize=HEAD, style='italic')
    yh = top - 1.6
    for name, x, w in cols:
        ax.text(x + w / 2 if name != 'PC' else x, yh, name, ha='center' if name != 'PC' else 'left', va='center', fontsize=HEAD, fontweight='bold')
    row_h = 0.85
    for i, r in enumerate(ho.itertuples()):
        y = yh - 1.0 - i * row_h
        ax.text(cols[0][1], y, f'Purkinje cell PC{r.PC}', ha='left', va='center', fontsize=BODY)
        for (name, x, w), s in zip(cols[1:5], (r.disc_sign_5W, r.disc_sign_18W, r.disc_sign_24W, r.disc_sign_30W)):
            ax.add_patch(plt.Rectangle((x + 0.05, y - 0.32), w - 0.1, 0.64, color=KI_COLOR if s == '+' else (WT_COLOR if s == '-' else '#bbbbbb'), alpha=0.85, lw=0))
            ax.text(x + w / 2, y, s.replace('-', MINUS), ha='center', va='center', color='white', fontsize=BODY, fontweight='bold')
        x, w = cols[5][1:]; ax.text(x + w / 2, y, r.predicted_sign.replace('-', MINUS), ha='center', va='center', fontsize=BODY, fontweight='bold')
        x, w = cols[6][1:]; ax.text(x + w / 2, y, f'{r.w12_diff_KI_minus_WT:+.3g}'.replace('-', MINUS), ha='center', va='center', fontsize=BODY)
        for (name, x, w), v in ((cols[7], r.sign_match), (cols[8], r.complete_separation_in_predicted_direction)):
            ax.text(x + w / 2, y, {'yes': '✓', 'no': '✗'}[v], ha='center', va='center', fontsize=MARK, color={'yes': YES, 'no': NO}[v])
        ax.plot([0.3, WIDTH_CM - 0.2], [y - row_h / 2, y - row_h / 2], color='#dddddd', lw=0.8)
    ax.text(0.3, 0.55, f'Orange + = KI higher, blue {MINUS} = KI lower. Predicted sign = majority of the four discovery signs.\n12 W difference = mean(KI) {MINUS} mean(WT) of three per-mouse means each.',
            ha='left', va='center', fontsize=NOTE, color='#444444', wrap=False)
    return fig


def main():
    ho, gate = load()
    fig = draw(ho); fig.canvas.draw()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = [t.get_fontsize() for t in texts]
    r = fig.canvas.get_renderer(); fw = fig.bbox.width
    overflow = [t.get_text() for t in texts if t.get_window_extent(r).x0 < 0 or t.get_window_extent(r).x1 > fw]
    assert not overflow, ('text runs past the figure edge', overflow)
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_S4_purkinje_only.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_S4_purkinje_only.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_S4_purkinje_only.svg', metadata={'Date': None}); plt.close(fig)
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    rec = dict(task='Supplementary S4 Purkinje-only', script_sha256=sha(Path(__file__)), inputs=dict(PINS), gate=gate,
               smallest_font_pt=min(sizes), font_sizes_used=sorted(set(sizes)), width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), dpi=300, matplotlib=matplotlib.__version__,
               outputs=outs, finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f: f.write(json.dumps(rec, indent=1) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS': f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(dict(gate=gate, smallest_font_pt=min(sizes), font_sizes=sorted(set(sizes)), width_cm=w_cm, height_cm=h_cm, **outs), indent=1))


if __name__ == '__main__':
    main()

#%% debug
