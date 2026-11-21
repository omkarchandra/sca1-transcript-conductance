
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Figure 2, panel A only: the channel heatmap (old Figure 2 panel B), rebuilt to the font rule (
# Plotting only. Source table: the S08 limma release table S08_S11_release_purkinje_bergmann.tsv.gz (40593a55), Purkinje rows. As in the accepted
# figure: colour = log2FC_decision (KI - WT), clipped at +-1.5; the number is printed where decision == HIT (BH q < 0.05); '*' = ambient_flagged_hit.
# The gene order, groups, gene/protein labels and ages are those of the accepted Figure_2 (final_figures_2.1_phase1_to_phase4/Figure_2.png 14 Oct).
# Colour scale: three-stop linear map #2b78d6 / #f0efeb / #e34948, read off the accepted figure's colour bar (mean RGB error 0.002).
#   layout: 17.5 cm wide, 13.5 cm tall, so the later forest plot (panel B) can sit below on the same page; built at that size, never rescaled.
#   font rule: ticks, cell numbers, gene labels 8 pt; group labels, axis/colour-bar labels 9 pt; panel letter 12 pt bold; note 8 pt; min font read back.
#   outputs: PNG 300 dpi, PDF, SVG (text kept as text), RECORD.json, SHA256SUMS into a NEW real_v01/.
# Run: python -I -B 06_figures/fig3_channel_heatmap.py
import hashlib, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/fig3_channel_heatmap'))
REL = ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz')
REL_SHA = '40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298'
AGES = (5, 18, 24, 30)
GROUPS = [('Pacemaking', [('Hcn1', 'HCN1')]),
          ('Spike initiation (Na⁺)', [('Scn8a', 'Nav1.6'), ('Scn1a', 'Nav1.1'), ('Scn2a', 'Nav1.2'), ('Scn1b', 'Navβ1'), ('Scn4b', 'Navβ4'), ('Fgf14', 'FGF14')]),
          ('Threshold (Kv1)', [('Kcna1', 'Kv1.1'), ('Kcna2', 'Kv1.2'), ('Kcnab1', 'Kvβ1')]),
          ('Low-threshold Ca²⁺, A-type K⁺', [('Cacna1g', 'Cav3.1'), ('Kcnd3', 'Kv4.3'), ('Kcnip1', 'KChIP1'), ('Dpp10', 'DPP10')]),
          ('Ca²⁺ entry at peak', [('Cacna1a', 'Cav2.1'), ('Cacna2d2', 'α2δ2')]),
          ('Repolarization', [('Kcnc3', 'Kv3.3'), ('Kcnma1', 'BK'), ('Kcnmb4', 'BK β4')]),
          ('After-hyperpolarization', [('Kcnn2', 'SK2'), ('Kcnq2', 'Kv7.2'), ('Kcnq3', 'Kv7.3'), ('Kcnq5', 'Kv7.5')]),
          ('Ca²⁺ store and extrusion', [('Itpr1', 'IP3R1'), ('Atp2a3', 'SERCA3'), ('Slc8a1', 'NCX1')])]
CMAP = LinearSegmentedColormap.from_list('fig2', ['#2b78d6', '#f0efeb', '#e34948']); CLIP = 1.5
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 13.5
TICK, LABEL, GROUP, NOTE, LETTER = 8, 9, 9, 8, 12
MINUS = '−'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    assert sha(REL) == REL_SHA, 'release table changed'
    t = pd.read_csv(REL, sep='\t', low_memory=False)
    genes = [g for _, gs in GROUPS for g, _ in gs]
    p = t[(t.population == 'purkinje') & t.dataset_symbol.isin(genes)]
    assert len(p) == len(genes) * len(AGES) and not p.duplicated(['dataset_symbol', 'age_weeks']).any()
    assert set(p.decision) <= {'HIT', 'TESTED_NOT_A_HIT'} and set(p.ambient_flagged_hit.astype(str)) <= {'True', 'False'}
    assert not (p.ambient_flagged_hit.astype(str).eq('True') & p.decision.ne('HIT')).any()
    key = p.set_index(['dataset_symbol', 'age_weeks'])
    val = np.array([[key.loc[(g, a), 'log2FC_decision'] for a in AGES] for g in genes], float)
    hit = np.array([[key.loc[(g, a), 'decision'] == 'HIT' for a in AGES] for g in genes])
    amb = np.array([[str(key.loc[(g, a), 'ambient_flagged_hit']) == 'True' for a in AGES] for g in genes])
    assert np.isfinite(val).all()
    return genes, val, hit, amb


def draw(genes, val, hit, amb):
    fig = plt.figure(figsize=(WIDTH_CM * CM, HEIGHT_CM * CM))
    # axes geometry in cm: heatmap columns start at x0, 1.75 cm per age; colour bar right of it
    x0, cw, top, bottom = 7.9, 1.75, 0.9, 1.75
    ax = fig.add_axes([x0 / WIDTH_CM, bottom / HEIGHT_CM, len(AGES) * cw / WIDTH_CM, (HEIGHT_CM - top - bottom) / HEIGHT_CM])
    n = len(genes); norm = Normalize(-CLIP, CLIP)
    for i in range(n):
        for j in range(len(AGES)):
            v = val[i, j]; ax.add_patch(plt.Rectangle((j + 0.04, i + 0.05), 0.92, 0.9, color=CMAP(norm(np.clip(v, -CLIP, CLIP))), lw=0))
            if hit[i, j]:
                s = f'{v:+.2f}'.replace('-', MINUS) + ('*' if amb[i, j] else '')
                ax.text(j + 0.5, i + 0.5, s, ha='center', va='center', fontsize=TICK, fontweight='bold', color='white' if abs(v) >= 0.75 else '#1b1b1b')
    ax.set_xlim(0, len(AGES)); ax.set_ylim(n, 0)
    ax.set_xticks(np.arange(len(AGES)) + 0.5, [f'{a} W' for a in AGES], fontsize=LABEL); ax.xaxis.tick_top(); ax.tick_params(axis='x', length=0, pad=3)
    labels = [f'{g} ({p})' for _, gs in GROUPS for g, p in gs]
    ax.set_yticks(np.arange(n) + 0.5, labels, fontsize=TICK); ax.tick_params(axis='y', length=0, pad=4)
    for s in ax.spines.values(): s.set_visible(False)
    # group brackets and labels, left of the gene labels
    i = 0; bx = -1.55  # in axes-data units (one age column = 1)
    for name, gs in GROUPS:
        lo, hi = i + 0.15, i + len(gs) - 0.15
        ax.plot([bx, bx], [lo, hi], color='#777777', lw=1.0, clip_on=False, transform=ax.transData)
        ax.text(bx - 0.1, (lo + hi) / 2, name, ha='right', va='center', fontsize=GROUP, color='#444444', clip_on=False)
        i += len(gs)
    cax = fig.add_axes([(x0 + len(AGES) * cw + 0.45) / WIDTH_CM, (bottom + 2.2) / HEIGHT_CM, 0.35 / WIDTH_CM, (HEIGHT_CM - top - bottom - 4.4) / HEIGHT_CM])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, ticks=[-1.5, -0.75, 0, 0.75, 1.5])
    cb.ax.set_yticklabels([f'{t:.2f}'.replace('-', MINUS) for t in (-1.5, -0.75, 0, 0.75, 1.5)], fontsize=TICK); cb.outline.set_visible(False)
    cb.set_label('KI ' + MINUS + ' WT log2 fold change', fontsize=LABEL)
    fig.text(0.35 / WIDTH_CM, 1 - 0.35 / HEIGHT_CM, 'A', fontsize=LETTER, fontweight='bold', va='top', ha='left')
    fig.text(0.35 / WIDTH_CM, 0.55 / HEIGHT_CM, 'SCA1 KI vs WT, Purkinje nuclei. Number shown = significant (BH q < 0.05); no number = not significant.\n'
             '* = significant but possibly ambient RNA. Colour clipped at ±1.5.', fontsize=NOTE, color='#444444', va='center', ha='left')
    return fig


def main():
    genes, val, hit, amb = load()
    fig = draw(genes, val, hit, amb); fig.canvas.draw()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = [t.get_fontsize() for t in texts]
    r = fig.canvas.get_renderer(); fw, fh = fig.bbox.width, fig.bbox.height
    overflow = [t.get_text() for t in texts if (e := t.get_window_extent(r)).x0 < -0.5 or e.x1 > fw + 0.5 or e.y0 < -0.5 or e.y1 > fh + 0.5]
    assert not overflow, ('text runs past the figure edge', overflow)
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_2A_channel_heatmap.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_2A_channel_heatmap.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_2A_channel_heatmap.svg', metadata={'Date': None}); plt.close(fig)
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    cells = {g: {f'{a}W': dict(log2FC_decision=float(val[i, j]), hit=bool(hit[i, j]), ambient=bool(amb[i, j])) for j, a in enumerate(AGES)} for i, g in enumerate(genes)}
    rec = dict(task='Figure 2 heatmap only (old panel B)', script_sha256=sha(Path(__file__)), release_table=dict(path=str(REL), sha256=REL_SHA),
               n_genes=len(genes), n_hit_cells=int(hit.sum()), n_ambient_cells=int(amb.sum()), colour_scale=dict(stops=['#2b78d6', '#f0efeb', '#e34948'], clip=CLIP), cells=cells,
               smallest_font_pt=min(sizes), font_sizes_used=sorted(set(sizes)), width_cm=round(w_cm, 3), height_cm=round(h_cm, 3), dpi=300, matplotlib=matplotlib.__version__,
               outputs=outs, finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f: f.write(json.dumps(rec, indent=1) + '\n')
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS': f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(dict(n_genes=len(genes), n_hit_cells=int(hit.sum()), smallest_font_pt=min(sizes), font_sizes=sorted(set(sizes)), width_cm=w_cm, height_cm=h_cm, **outs), indent=1))


if __name__ == '__main__':
    main()

#%% debug
