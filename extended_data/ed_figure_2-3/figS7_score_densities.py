#!/usr/bin/env python3
"""Extended Data Figure 2-3: per-mouse Purkinje lane-score densities, from the E2 occupancy inputs."""
from __future__ import annotations

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
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.stats import gaussian_kde

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/ed_figure_2-3'))
AGES = (5, 18, 24, 30)
W12_DIRECTION_AGE = 18
P90 = 0.90
WIDTH_CM, HEIGHT_CM = 17.5, 6.5
LANES = {5: (5, 6, 11), 18: (5, 6, 11, 30), 24: (5, 6, 11, 30), 30: (5, 6, 30)}
INPUTS = {
    'nucleus_pcs': ('extended_data/ed_figure_2-3/inputs/nucleus_PCs.npy', 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b'),
    'nucleus_metadata': ('extended_data/ed_figure_2-3/inputs/nucleus_metadata.tsv.gz', '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721'),
    'unit_metadata': ('data/snrnaseq/unit_metadata.tsv', 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657'),
    'e2_source': ('01_snrnaseq/occupancy/e2_occupancy_lomo.py', '7531de39c4b65d4ec77f9359cf1eeef99634c9f809805915a708a45cc53b5e95'),
    'e2_record': ('extended_data/ed_figure_2-3/inputs/E2_RECORD.json', '05fd17596643e5cea19b9aaec5d324101dd0953e9fe3ad6bc55df6679bef4e7a'),
    'e2a': ('extended_data/ed_figure_2-3/inputs/E2a_LOMO.tsv', '511b2559c5c1c99269f4ca6b01a2febfa5ad6b2ac98ed55de195a78124f195b2'),
    'e2b': ('extended_data/ed_figure_2-3/inputs/E2b_W12_SHARES.tsv', 'e0b98e9e625bd9e4fdd86ce32bcb9bf2722e90288a9293dac41adcf52c5cacdd'),
    'far_occupancy': ('results/far_set_ratios/real_v02/OCCUPANCY.tsv', '2823a594148d1b625a4938d3cd8ecb12112a86478fd0b0535b1cf2bcaebf9b1e'),
    'w12_pcs': ('extended_data/ed_figure_2-3/inputs/Purkinje_cell_PCs_12w.npy', '3c2eb711b5f3c6777c2f7deb43789337633132194bb09d168544d312ddcfb310'),
    'w12_rows': ('extended_data/ed_figure_2-3/inputs/Purkinje_cell_projection_12w.tsv.gz', '26161c549f37faa58567eca24064603aedec955c5d875996754506e1359d15e4'),
}
WT_COLORS = ('#4D4D4D', '#858585', '#B5B5B5')
KI_COLORS = ('#9E2A2B', '#D1495B', '#F28E63')


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def token_tol(token: str, a: float, b: float) -> float:
    half_last = float(abs(Decimal(1).scaleb(Decimal(token).as_tuple().exponent)) / 2)
    return max(1e-9 * max(1.0, abs(a), abs(b)), half_last)


def read_table(key: str, dtype=None) -> pd.DataFrame:
    return pd.read_csv(ROOT / INPUTS[key][0], sep='\t', comment='#', dtype=dtype)


def main() -> None:
    # A new output folder; all generated products must be absent.
    OUT.mkdir(parents=True, exist_ok=False)
    for name in ('SCORES.tsv', 'GATE.tsv', 'RECORD.json', 'Figure_S7.png', 'Figure_S7.svg', 'Figure_S7.pdf', 'SHA256SUMS'):
        assert not (OUT / name).exists(), f'exclusive-create product exists: {name}'
    for key, (rel, expected) in INPUTS.items():
        actual = sha(ROOT / rel)
        assert actual == expected, (key, actual, expected)

    X = np.load(ROOT / INPUTS['nucleus_pcs'][0])
    nm = read_table('nucleus_metadata')
    um = read_table('unit_metadata', dtype={'animal_id': str})
    assert X.ndim == 2 and X.shape[0] == len(nm) and X.shape[1] >= 30 and np.isfinite(X).all()
    assert nm.obs_name.is_unique and (um.index.to_numpy() == um.unit_index.to_numpy()).all()
    pk = um[um.cell_type_raw == 'Purkinje cell'].copy()
    assert len(pk) == 24 and set(pk.age_weeks) == set(AGES) and pk.groupby(['age_weeks', 'genotype']).size().eq(3).all()
    cell = nm.merge(pk[['unit_index', 'animal_id', 'age_weeks', 'genotype']], on='unit_index', how='left', validate='many_to_one')
    assert cell.animal_id.notna().all() and (cell.obs_name.to_numpy() == nm.obs_name.to_numpy()).all()
    cell.animal_id = cell.animal_id.astype(str)

    environment = {'interpreter': sys.executable, 'scipy': scipy.__version__, 'numpy': np.__version__, 'matplotlib': matplotlib.__version__}
    e2_record = json.loads((ROOT / INPUTS['e2_record'][0]).read_text())
    directions = {}
    for age in AGES:
        vector = e2_record['full_unit_vector'][str(age)]
        u = np.array([vector[f'PC{k}'] for k in LANES[age]], dtype=float)
        assert abs(np.linalg.norm(u) - 1.0) < 1e-12
        directions[age] = u

    e2 = read_table('e2a', dtype=str)
    far = read_table('far_occupancy', dtype=str)
    e2b = read_table('e2b', dtype=str)
    gate_rows, score_parts = [], []
    discovery_raw, wt_norm = {}, {}

    # Exact scoring arithmetic copied from accepted E2 source lines 165-166:
    # X[rows][:, [k - 1 for k in EXPECTED_LANES[a]]] @ u
    for age in AGES:
        age_mask = cell.age_weeks.astype(int).eq(age).to_numpy()
        raw = X[age_mask][:, [k - 1 for k in LANES[age]]] @ directions[age]
        age_cell = cell.loc[age_mask, ['obs_name', 'animal_id', 'genotype']].reset_index(drop=True)
        assert len(raw) == len(age_cell)
        wt_mask = age_cell.genotype.eq('WT').to_numpy()
        wt_raw = raw[wt_mask]
        mu, sd = float(wt_raw.mean()), float(wt_raw.std(ddof=0))
        assert sd > 0
        threshold = float(np.quantile(wt_raw, P90))
        discovery_raw[age] = (age_cell, raw, threshold)
        wt_norm[age] = (mu, sd)
        z = (raw - mu) / sd
        score_parts.append(pd.DataFrame({
            'nucleus_id': age_cell.obs_name, 'mouse': age_cell.animal_id,
            'genotype': age_cell.genotype, 'age': age, 'direction_age': age,
            'held_out': False, 'lane_score': z, 'lane_score_raw': raw,
        }))
        for mouse in age_cell.animal_id.drop_duplicates():
            mask = age_cell.animal_id.eq(mouse).to_numpy()
            n = int(mask.sum()); count = int((raw[mask] > threshold).sum()); share = count / n
            ae = e2[(e2.age.astype(int) == age) & (e2.mouse == mouse)].iloc[0]
            af = far[(far.age.astype(int) == age) & (far.mouse == mouse)].iloc[0]
            e_thr, e_count, e_share = float(ae.full_p90), int(ae.full_n_beyond), float(ae.full_share)
            f_count, f_share = int(af.n_beyond_WT_p90), float(af.share_beyond_WT_p90)
            thr_tol = token_tol(ae.full_p90, threshold, e_thr)
            e_share_tol = token_tol(ae.full_share, share, e_share)
            f_share_tol = token_tol(af.share_beyond_WT_p90, share, f_share)
            ok = (count == e_count == f_count and abs(threshold-e_thr) <= thr_tol
                  and abs(share-e_share) <= e_share_tol and abs(share-f_share) <= f_share_tol)
            gate_rows.append(dict(panel_age=age, direction_age=age, held_out=False, mouse=mouse,
                                  genotype=age_cell.loc[mask, 'genotype'].iloc[0], n_nuclei=n,
                                  recomputed_p90=threshold, accepted_e2_p90=e_thr, p90_tolerance=thr_tol,
                                  recomputed_n_beyond=count, accepted_e2_n_beyond=e_count, accepted_far_n_beyond=f_count,
                                  recomputed_share=share, accepted_e2_share=e_share, e2_share_tolerance=e_share_tol,
                                  accepted_far_share=f_share, far_share_tolerance=f_share_tol, pass_gate=ok))

    # Held-out 12-week nuclei on the accepted 18-week full-data direction.
    X12 = np.load(ROOT / INPUTS['w12_pcs'][0])
    r12 = read_table('w12_rows', dtype={'animal_id': str})
    assert X12.shape == (len(r12), 30) and np.isfinite(X12).all()
    raw12 = X12[:, [k - 1 for k in LANES[W12_DIRECTION_AGE]]] @ directions[W12_DIRECTION_AGE]
    mu18, sd18 = wt_norm[W12_DIRECTION_AGE]
    threshold18 = discovery_raw[W12_DIRECTION_AGE][2]
    score_parts.append(pd.DataFrame({
        'nucleus_id': r12.obs_name, 'mouse': r12.animal_id.astype(str), 'genotype': r12.genotype,
        'age': 12, 'direction_age': W12_DIRECTION_AGE, 'held_out': True,
        'lane_score': (raw12-mu18)/sd18, 'lane_score_raw': raw12,
    }))
    for mouse in r12.animal_id.astype(str).drop_duplicates():
        mask = r12.animal_id.astype(str).eq(mouse).to_numpy()
        n = int(mask.sum()); count = int((raw12[mask] > threshold18).sum()); share = count/n
        ae = e2b[(e2b.reference_age.astype(int) == W12_DIRECTION_AGE) & (e2b.mouse == mouse)].iloc[0]
        e_thr, e_count, e_share = float(ae.ref_p90), int(ae.n_beyond), float(ae.share)
        thr_tol = token_tol(ae.ref_p90, threshold18, e_thr)
        share_tol = token_tol(ae.share, share, e_share)
        ok = count == e_count and abs(threshold18-e_thr) <= thr_tol and abs(share-e_share) <= share_tol
        gate_rows.append(dict(panel_age=12, direction_age=W12_DIRECTION_AGE, held_out=True, mouse=mouse,
                              genotype=r12.loc[mask, 'genotype'].iloc[0], n_nuclei=n,
                              recomputed_p90=threshold18, accepted_e2_p90=e_thr, p90_tolerance=thr_tol,
                              recomputed_n_beyond=count, accepted_e2_n_beyond=e_count, accepted_far_n_beyond=np.nan,
                              recomputed_share=share, accepted_e2_share=e_share, e2_share_tolerance=share_tol,
                              accepted_far_share=np.nan, far_share_tolerance=np.nan, pass_gate=ok))

    gate = pd.DataFrame(gate_rows)
    assert len(gate) == 30 and gate.mouse.nunique() == 30
    gate.to_csv(OUT / 'GATE.tsv', sep='\t', index=False, float_format='%.12g')
    if not gate.pass_gate.all():
        failed = gate.loc[~gate.pass_gate].to_dict(orient='records')
        rec = {'status': 'GATE_FAILED_NO_FIGURE', 'failed': failed, 'inputs': {rel: h for rel, h in INPUTS.values()},
               'script_sha256': sha(Path(__file__)), 'environment': environment, 'finished_utc': datetime.now(timezone.utc).isoformat(timespec='seconds')}
        (OUT/'RECORD.json').write_text(json.dumps(rec, indent=2)+'\n')
        raise SystemExit('gate failed; no score table or figure created')

    scores = pd.concat(score_parts, ignore_index=True)
    scores = scores.sort_values(['age', 'genotype', 'mouse', 'nucleus_id']).reset_index(drop=True)
    scores.to_csv(OUT / 'SCORES.tsv', sep='\t', index=False, float_format='%.12g')

    matplotlib.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'axes.labelsize':9,
                                'axes.titlesize':9, 'xtick.labelsize':8, 'ytick.labelsize':8,
                                'legend.fontsize':8, 'pdf.fonttype':42, 'svg.fonttype':'none',
                                'axes.linewidth':0.7})
    panels = [(5, False), (12, True), (18, False), (24, False), (30, False)]
    fig, axes = plt.subplots(1, 5, figsize=(WIDTH_CM/2.54, HEIGHT_CM/2.54), sharey=False)
    fig.subplots_adjust(left=.075, right=.975, top=.60, bottom=.22, wspace=.28)
    global_min = float(scores.lane_score.quantile(.001)); global_max = float(scores.lane_score.quantile(.999))
    pad = .14*(global_max-global_min); xmin, xmax = global_min-pad, global_max+pad
    grid = np.linspace(xmin, xmax, 500)
    legend_rows = []
    letter_artists = []
    for i, (age, held) in enumerate(panels):
        ax = axes[i]
        panel = scores[(scores.age == age) & (scores.held_out == held)]
        gpanel = gate[(gate.panel_age == age) & (gate.held_out == held)]
        direction_age = W12_DIRECTION_AGE if held else age
        mu, sd = wt_norm[direction_age]
        threshold_z = (discovery_raw[direction_age][2]-mu)/sd
        wt_mice = sorted(panel.loc[panel.genotype=='WT','mouse'].unique())
        ki_mice = sorted(panel.loc[panel.genotype=='KI','mouse'].unique())
        assert len(wt_mice)==len(ki_mice)==3
        ymax = 0.0
        for mice, colors, genotype in ((wt_mice, WT_COLORS, 'WT'), (ki_mice, KI_COLORS, 'KI')):
            for mouse, color in zip(mice, colors):
                values = panel.loc[panel.mouse == mouse, 'lane_score'].to_numpy(float)
                density = gaussian_kde(values)(grid)
                ymax=max(ymax,float(density.max()))
                ax.plot(grid,density,color=color,lw=1.15 if genotype=='KI' else 1.0,alpha=.98)
        ax.axvline(threshold_z,color='#222222',ls='--',lw=.8)
        for genotype, mice, colors in (('WT', wt_mice, WT_COLORS), ('KI', ki_mice, KI_COLORS)):
            for mouse, color in zip(mice, colors):
                share=float(gpanel.loc[gpanel.mouse==mouse,'recomputed_share'].iloc[0])
                legend_rows.append(dict(panel_age=age, held_out=held, direction_age=direction_age, genotype=genotype, mouse=mouse, occupancy=share, color=color))
        for j,(mouse,color) in enumerate(zip(ki_mice,KI_COLORS)):
            share=float(gpanel.loc[gpanel.mouse==mouse,'recomputed_share'].iloc[0])
            ax.text(.97,.94-j*.14,f'{share:.0%}',transform=ax.transAxes,ha='right',va='center',fontsize=8,color=color)
        ax.set_xlim(xmin,xmax); ax.set_ylim(0,ymax*1.12)
        ax.set_xticks([t for t in np.arange(np.ceil(xmin/5)*5, xmax, 5) if xmin + .08*(xmax-xmin) < t < xmax - .08*(xmax-xmin)])
        ax.set_title('12 weeks (held-out)\n18-week axis direction' if held else f'{age} weeks',pad=3,fontsize=9)
        letter_artists.append(ax.text(-.14,1.40,chr(65+i),transform=ax.transAxes,fontsize=12,fontweight='bold',ha='left',va='bottom'))
        ax.spines[['top','right']].set_visible(False); ax.tick_params(width=.7,length=3)
        ax.grid(axis='y',color='#E7E9ED',lw=.6); ax.set_axisbelow(True)
        if i==0: ax.set_ylabel('Density')
        else: ax.set_yticklabels([]); ax.set_ylabel('')
    fig.supxlabel('Score along the axis direction (WT-standardised units)',fontsize=9,y=.055)
    legend_table = pd.DataFrame(legend_rows)
    legend_table.to_csv(OUT/'LEGEND.tsv',sep='\t',index=False,float_format='%.12g')
    handles=[Line2D([],[],color=WT_COLORS[0],lw=1.2,label='WT mice'),
             Line2D([],[],color=KI_COLORS[1],lw=1.3,label='KI mice (label = shifted fraction, beyond the WT 90th percentile)'),
             Line2D([],[],color='#222222',lw=.8,ls='--',label='WT 90th percentile')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.995),ncol=3,frameon=False,fontsize=8)
    fig.canvas.draw()
    title_artists = [ax.title for ax in axes]
    letter_title_overlap = [(letter.get_text(), title.get_text()) for letter, title in zip(letter_artists, title_artists) if letter.get_window_extent(fig.canvas.get_renderer()).overlaps(title.get_window_extent(fig.canvas.get_renderer()))]
    assert not letter_title_overlap, ('panel letter/title overlap', letter_title_overlap)
    texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    min_font=min(float(t.get_fontsize()) for t in texts); assert min_font>=8
    r=fig.canvas.get_renderer(); fw,fh=fig.bbox.width,fig.bbox.height
    overflow=[t.get_text() for t in texts if (b:=t.get_window_extent(r)).x0<-.5 or b.x1>fw+.5 or b.y0<-.5 or b.y1>fh+.5]
    assert not overflow, ('text overflow',overflow)
    for ext,meta in [('png',{'Software':None}),('svg',{'Date':None}),('pdf',{'CreationDate':None,'ModDate':None})]:
        fig.savefig(OUT/f'Figure_S7.{ext}',dpi=300,metadata=meta)
    plt.close(fig)

    generated=[OUT/'GATE.tsv',OUT/'SCORES.tsv',OUT/'LEGEND.tsv',OUT/'Figure_S7.png',OUT/'Figure_S7.svg',OUT/'Figure_S7.pdf']
    record={
        'status':'PASS','task':'Extended Data Figure 2-3 (drawn as Supplementary Figure S7): per-mouse score densities',
        'script_sha256':sha(Path(__file__)),
        'environment':environment,
        'accepted_scoring_source':'exact arithmetic of 01_snrnaseq/occupancy/e2_occupancy_lomo.py (lines 165-166 as run)',
        'inputs':{rel:h for rel,h in INPUTS.values()},
        'full_unit_vectors_read_from_E2_RECORD':{str(a):{f'PC{k}':float(v) for k,v in zip(LANES[a],directions[a])} for a in AGES},
        'lomo_direction_columns_used':False,
        'gate':{'rows':len(gate),'passed':int(gate.pass_gate.sum()),'failed':int((~gate.pass_gate).sum()),
                'rule':'abs(a-b) <= max(1e-9*max(1,abs(a),abs(b)), half unit in last stored digit); counts exact'},
        'scores':{'rows':len(scores),'discovery_rows':int((~scores.held_out).sum()),'heldout_12w_rows':int(scores.held_out.sum()),
                  'lane_score':'raw full-direction score standardised by pooled discovery WT mean and population SD (ddof=0) within direction age'},
        'figure':{'panels':5,'width_cm':WIDTH_CM,'height_cm':HEIGHT_CM,'dpi':300,'minimum_font_pt':min_font,
                  'svg_text_as_text':matplotlib.rcParams['svg.fonttype']=='none','visual_overflow_check':overflow},
        'outputs':{p.name:sha(p) for p in generated},
        'finished_utc':datetime.now(timezone.utc).isoformat(timespec='seconds'),
    }
    (OUT/'RECORD.json').write_text(json.dumps(record,indent=2)+'\n')
    all_files=[Path(__file__).resolve()]+generated+[OUT/'RECORD.json']
    with (OUT/'SHA256SUMS').open('x') as f:
        for p in all_files: f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
