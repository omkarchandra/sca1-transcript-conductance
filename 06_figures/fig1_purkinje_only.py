
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Figure 1, Purkinje-only.. Plotting only, from accepted tables; no new numbers
# beyond the plotted mean of three per-mouse values (dashed line), as in the original Figure 1 panel A (the Phase 4 figure script figure_lane_means).
#   values: per-mouse discovery PC means at 5/18/24/30 W from the frozen S12 T3 mouse table; the six 12 W mice (out-of-sample projection, no refit)
#           from the accepted 12 W PC-means table, cross-checked against HELDOUT_BY_PC.tsv (f32ed294) w12 means to its printed precision.
#   gate:   Purkinje sign matches 4 of 4 and complete separation 3 of 4 (PC5, PC6, PC30) in HELDOUT_BY_PC.tsv, else stop.
#   font rule: built at final width 17.5 cm (2 x 2), never rescaled (no bbox_inches='tight'); ticks 8 pt, axis labels and titles 9 pt, legend 8 pt,
#           panel letters 12 pt bold; the smallest font size is read back from the figure object and written to RECORD.json.
# v02: v01 plus an SVG save with text kept as text (svg.fonttype = 'none'); output folder real_v02/.
# Run: python -I -B 06_figures/fig1_purkinje_only.py
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/fig1_purkinje_only'))
PCS = (5, 6, 11, 30); SIX = ('2061', '2062', '22018', '2063', '2065', '22019'); AGES = (5, 18, 24, 30)
T3 = _release_paths.expand('${REPO_ROOT}/results/s12_lanes/mouse_table.tsv'); M12 = _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Purkinje_cell_PC_means_by_mouse_12w.tsv')
HO = _release_paths.expand('${REPO_ROOT}/results/heldout/HELDOUT_BY_PC.tsv'); MAN = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv')
PINS = {T3: '94cdfb6f7213dfee8d9ad01f30c0b6f99404335007a234cb2c4dc7b6d3938b00', M12: '4613c11497fe055297a540595a7b855849ba9595e0ca78deabcb15a4c7ea9158',
        HO: 'f32ed29456c1dbafae14369dcba9e358a90ae0d240b3d94a819f4e377ed73d5b', MAN: '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'}
WT_COLOR, KI_COLOR = '#0072B2', '#D55E00'
CM = 1 / 2.54; WIDTH_CM, HEIGHT_CM = 17.5, 14.0
TICK, LABEL, TITLE, LEGEND, LETTER = 8, 9, 9, 8, 12


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def load():
    for rel, h in PINS.items():
        assert sha(ROOT / rel) == h, ('pinned file changed', rel)
    man = pd.read_csv(ROOT / MAN, sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    t3 = pd.read_csv(ROOT / T3, sep='\t', dtype={'animal_id': str}); d = t3[t3.population == 'Purkinje_cell']
    assert len(d) == 24 and sorted(d.age_weeks.unique()) == list(AGES) and (d.genotype.to_numpy() == man.genotype.loc[d.animal_id].to_numpy()).all()
    assert (d.groupby(['age_weeks', 'genotype']).size() == 3).all()
    m12 = pd.read_csv(ROOT / M12, sep='\t', dtype={'animal_id': str}).set_index('animal_id'); assert sorted(m12.index) == sorted(SIX)
    assert list(m12.genotype.loc[list(SIX)]) == [man.genotype[a] for a in SIX]
    ho = pd.read_csv(ROOT / HO, sep='\t'); ho = ho[ho.population == 'Purkinje_cell'].set_index('PC'); assert tuple(ho.index) == PCS
    n_match, sep = int((ho.sign_match == 'yes').sum()), tuple(int(k) for k in ho.index[ho.complete_separation_in_predicted_direction == 'yes'])
    if n_match != 4 or sep != (5, 6, 30):
        raise SystemExit(f'STOP: held-out gate failed: sign matches {n_match}/4, complete separation {sep}')
    vals, dev = {}, 0.0
    for k in PCS:
        for col in ('w12_WT_means', 'w12_KI_means'):
            for item in ho.loc[k, col].split(';'):
                a, v = item.split(':'); dev = max(dev, abs(float(v) - m12.loc[a, f'PC{k}']) / max(1.0, abs(float(v))))
        a_ = pd.DataFrame(dict(animal_id=d.animal_id, genotype=d.genotype, age_weeks=d.age_weeks, value=d[f'discovery_PC{k}_mean'], sample_set='in-sample'))
        b_ = pd.DataFrame(dict(animal_id=list(SIX), genotype=[man.genotype[x] for x in SIX], age_weeks=12, value=[m12.loc[x, f'PC{k}'] for x in SIX], sample_set='out-of-sample'))
        vals[k] = pd.concat([a_, b_], ignore_index=True)
    assert dev < 1e-5, ('12 W table disagrees with HELDOUT_BY_PC.tsv', dev)
    return vals, ho, dict(sign_matches=n_match, complete_separation=list(sep), max_rel_dev_12w_vs_heldout=dev)


def draw(vals, ho):
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_CM * CM, HEIGHT_CM * CM), layout='constrained')
    for ax, k, letter in zip(axes.ravel(), PCS, 'ABCD'):
        d = vals[k]; r = ho.loc[k]
        ax.axvspan(10.5, 13.5, color='#dddddd', alpha=0.6, lw=0)
        for geno, color, dx in (('WT', WT_COLOR, -0.45), ('KI', KI_COLOR, 0.45)):
            g = d[d.genotype == geno]; means = g.groupby('age_weeks').value.mean()
            ax.plot(means.index, means.values, ls='--', lw=1.2, color=color, alpha=0.8, zorder=1)
            s = g[g.sample_set == 'in-sample']; ax.scatter(s.age_weeks + dx, s.value, s=28, color=color, edgecolor='white', lw=0.5, zorder=3)
            o = g[g.sample_set == 'out-of-sample']; ax.scatter(o.age_weeks + dx, o.value, s=40, marker='D', facecolor='white', edgecolor=color, lw=1.4, zorder=4)
        ax.set_xticks([5, 12, 18, 24, 30]); ax.tick_params(labelsize=TICK); ax.set_xlabel('age (weeks)', fontsize=LABEL)
        ax.set_ylabel('mean PC score per mouse', fontsize=LABEL)
        ax.set_title(f"Purkinje cell PC{k}\npred. {r.predicted_sign.replace('-', '−')} · match {r.sign_match} · sep. {r.complete_separation_in_predicted_direction}", fontsize=TITLE)
        ax.text(-0.22, 1.13, letter, transform=ax.transAxes, fontsize=LETTER, fontweight='bold', va='bottom', ha='left')
    handles = [plt.Line2D([], [], marker='o', ls='', color=WT_COLOR, label='WT, 24-mouse set (in-sample)'), plt.Line2D([], [], marker='o', ls='', color=KI_COLOR, label='KI, 24-mouse set (in-sample)'),
               plt.Line2D([], [], marker='D', ls='', mfc='white', mec=WT_COLOR, mew=1.4, label='WT, 12 W held out (projected, no refit)'), plt.Line2D([], [], marker='D', ls='', mfc='white', mec=KI_COLOR, mew=1.4, label='KI, 12 W held out (projected, no refit)'),
               plt.Line2D([], [], ls='--', color='grey', label='dashed: mean of the 3 mice per genotype'), plt.Rectangle((0, 0), 1, 1, color='#dddddd', alpha=0.6, label='grey band: 12 W, not used to find the components')]
    fig.legend(handles=handles, loc='outside lower center', ncol=2, fontsize=LEGEND, frameon=False,
               title='pred. = predicted KI − WT sign (majority over 5/18/24/30 W); match = 12 W sign matched the prediction;\nsep. = complete 3 v 3 separation of the 12 W mice in the predicted direction',
               title_fontsize=LEGEND)
    return fig


def main():
    vals, ho, gate = load()
    fig = draw(vals, ho); fig.canvas.draw()
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    w_cm, h_cm = (x * 2.54 for x in fig.get_size_inches())
    OUT.mkdir(parents=False, exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42; matplotlib.rcParams['svg.fonttype'] = 'none'
    fig.savefig(OUT / 'Figure_1_purkinje_only.png', dpi=300, metadata={'Software': None}); fig.savefig(OUT / 'Figure_1_purkinje_only.pdf', metadata={'CreationDate': None})
    fig.savefig(OUT / 'Figure_1_purkinje_only.svg', metadata={'Date': None}); plt.close(fig)
    outs = {p.name: sha(p) for p in sorted(OUT.iterdir())}
    rec = dict(task='Figure 1 Purkinje-only', script_sha256=sha(Path(__file__)), inputs={rel: h for rel, h in PINS.items()},
               lane_means_table=dict(path=T3, sha256=PINS[T3]), heldout_table=dict(path=HO, sha256=PINS[HO]), gate=gate,
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
