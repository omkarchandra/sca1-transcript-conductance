"""Figure 4: accepted replication estimates only; no refitting.

Called by fig4_replication.py for the standalone replication figure.
Six strips preserve legible 9 pt labels at the requested 17.5 cm width.
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
BASE = ROOT / _release_paths.expand('${REPO_ROOT}/results')
WIDTH, HEIGHT = 17.5, 17.0
GENES = ['Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g',
         'Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2']
MANIFEST = ROOT / _release_paths.expand('${REPO_ROOT}/data/INPUT_MANIFEST.tsv')
SOURCES = {
    'a1': (BASE / _release_paths.expand('${REPO_ROOT}/results/a1/primary/NINE_GENES.tsv'),
           'ed5f07e9268f65bc326df9b94cdbdf41fa82d9a1df2e28bc0041a439aa7a49ff'),
    'a1_reading': (BASE / _release_paths.expand('${REPO_ROOT}/results/a1/primary/READING.tsv'),
                   'd7d61c71629ebc14724a2eb08778d81da4de0f6dfcf577ea77675e9cce34ec0f'),
    'expl': (BASE / _release_paths.expand('${REPO_ROOT}/results/a1/exploratory/NINE_GENES.tsv'),
             'dfe88a877583165929bfcbe8d04b1bc075ef1bed867dcd89f4c20f5feb8235c2'),
    'expl_reading': (BASE / _release_paths.expand('${REPO_ROOT}/results/a1/exploratory/READING.tsv'),
                     '399bc2bcc0c4d98b00630e0224a0705aa9baa213dfc3c28b19aa0efdbc77eb44'),
    'p2': (BASE / _release_paths.expand('${REPO_ROOT}/results/p2/NINE_GENES.tsv'),
           '22cba0dabd7e08b38b47cc4393ea77fbbc7aab96fd45fb375a00ffb8533cfd55'),
    'e3': (Path(_release_paths.expand('${REPO_ROOT}/data/snrnaseq/E3_LIMMA_RESULTS.tsv')),
           'eab5ccbedff4c8b340ed6b6a9c43920bcfe7aa308387a2702f3cd8ab4c101e67'),
}
DATASETS = [
    ('Discovery Purkinje nuclei, 12 weeks (held-out mice)', 'replicates (5/5, 9/9)', 'e3', None, 12),
    ('ATXN1[82Q] Purkinje nuclei, 12 weeks', 'replicates (5/5, 9/9)', 'a1', 'GSE215336', 12),
    ('SCA7 Purkinje nuclei, 8.5 weeks', 'exploratory (9/9)', 'expl', 'GSE269430_SCA7', 8.5),
    ('SCA7 Purkinje nuclei, 5 weeks', 'exploratory (9/9)', 'expl', 'GSE269430_SCA7', 5),
    ('Knock-in bulk cerebellum, 12 weeks', 'not met (3/5)', 'p2', 'GSE122099', 12),
    ('Human SCA1 Purkinje nuclei, 13 donors (primary)', 'inconclusive (2/5)', 'a1', 'GSE246183_human', 0),
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load():
    manifest_text = MANIFEST.read_text()
    tables = {}
    for key, (path, expected) in SOURCES.items():
        assert expected in manifest_text, (key, 'not pinned in MANIFEST_A')
        assert sha(path) == expected, (key, 'source changed')
        tables[key] = pd.read_csv(path, sep='\t', comment='#', low_memory=False)
    records = []
    for title, verdict, source, dataset, age in DATASETS:
        t = tables[source]
        if source == 'e3':
            selected = t[(t.population == 'Purkinje cell') & (t.age_weeks == age)]
            gene, value, lo, hi = 'dataset_symbol', 'log2FC', 'CI_L', 'CI_R'
        elif source == 'p2':
            selected = t[(t.dataset == dataset) & (t.age == age) & (t.role == 'PRIMARY') & (t.kcnc3 == 'merged')]
            gene, value, lo, hi = 'gene', 'log2fc_mut_vs_wt', 'ci_lo', 'ci_hi'
        else:
            selected = t[(t.dataset == dataset) & (t.age == age)]
            gene, value, lo, hi = 'gene', 'logFC', 'CI.L', 'CI.R'
            reading = tables['expl_reading' if source == 'expl' else 'a1_reading']
            r = reading[(reading.dataset == dataset) & (reading.age == age)]
            assert len(r) == 1
        selected = selected[selected[gene].isin(GENES)]
        assert len(selected) == 9 and not selected[gene].duplicated().any()
        selected = selected.set_index(gene).loc[GENES]
        vals = selected[[value, lo, hi]].to_numpy(float)
        assert np.isfinite(vals).all()
        assert (vals[:, 1] <= vals[:, 0]).all() and (vals[:, 0] <= vals[:, 2]).all()
        negatives = (int((vals[:5, 0] < 0).sum()), int((vals[:, 0] < 0).sum()))
        assert negatives == ((3, 6) if source == 'p2' else (2, 5) if dataset == 'GSE246183_human' else (5, 9))
        if source not in ('e3', 'p2'):
            assert negatives == (int(r.iloc[0].n_channel_negative), int(r.iloc[0].n_nine_negative))
        for i, g in enumerate(GENES):
            row = selected.loc[g]
            records.append(dict(dataset=title, verdict=verdict, source=source, gene=g,
                                log2FC=float(vals[i, 0]), ci_low=float(vals[i, 1]), ci_high=float(vals[i, 2]),
                                source_row_symbol=str(row.get('row_symbol', g)),
                                ci_method=str(row.get('ci_method', 'limma 95% CI'))))
    return pd.DataFrame(records)


def draw(fig, table, total_height=HEIGHT):
    # All geometry is in physical cm. Panel B occupies the bottom HEIGHT cm.
    fig.text(.35 / WIDTH, (HEIGHT - .35) / total_height, 'B', fontsize=12, fontweight='bold', va='top')
    fig.text(1.3 / WIDTH, (HEIGHT - .4) / total_height,
             'Replication: KI − WT or disease − control log2 fold change', fontsize=9, va='top')
    fig.text(1.3 / WIDTH, (HEIGHT - .95) / total_height,
             'Five channel genes', fontsize=9, color='#286b9d', va='top')
    fig.text(10.1 / WIDTH, (HEIGHT - .95) / total_height,
             'Four other genes', fontsize=9, color='#765483', va='top')
    artists = []
    for i, (title, verdict, source, dataset, age) in enumerate(DATASETS):
        top = HEIGHT - 1.95 - i * 2.28
        bottom, ax_height = top - 1.35, 1.35
        ax = fig.add_axes([1.3 / WIDTH, bottom / total_height, 15.65 / WIDTH, ax_height / total_height])
        rows = table[table.dataset == title].set_index('gene').loc[GENES]
        v, lo, hi = (rows[c].to_numpy(float) for c in ('log2FC', 'ci_low', 'ci_high'))
        ax.axvspan(-.5, 4.5, color='#286b9d', alpha=.055, lw=0)
        ax.axvspan(4.5, 8.5, color='#765483', alpha=.055, lw=0)
        ax.axvline(4.5, color='#aaaaaa', lw=.65)
        ax.axhline(0, color='#333333', lw=.8, zorder=2)
        bars = ax.bar(np.arange(9), v, width=.62, color=['#347eae'] * 5 + ['#91739b'] * 4,
                      yerr=np.vstack([v - lo, hi - v]), capsize=2,
                      error_kw={'elinewidth': .8, 'capthick': .8, 'ecolor': '#333333'}, zorder=3)
        artists.extend(zip(bars.patches, rows.to_dict('records')))
        ax.set_xlim(-.5, 8.5)
        # Shared x axis; identical gene positions in all six strips.
        ax.set_xticks(np.arange(9), [g + ('†' if source == 'p2' and g == 'Kcnc3' else '') for g in GENES], fontsize=9)
        ax.tick_params(axis='x', length=0, pad=3, labelbottom=i == 5)
        ax.tick_params(axis='y', labelsize=8, length=2, pad=2)
        ax.set_ylim(-2.2, 3.2)
        ax.set_yticks([-2, 0, 2])
        ax.spines[['top', 'right', 'bottom']].set_visible(False)
        ax.spines['left'].set_color('#aaaaaa')
        fig.text(1.3 / WIDTH, (top + .22) / total_height, title, fontsize=9, va='bottom')
        fig.text(16.95 / WIDTH, (top + .22) / total_height, verdict, fontsize=9, ha='right', va='bottom', fontweight='bold')
        assert lo.min() >= -2.2 and hi.max() <= 3.2, (title, 'error bar clipped')
        # Check the errorbar segment coordinates actually handed to Matplotlib.
        segments = bars.errorbar.lines[2][0].get_segments()
        for j, seg in enumerate(segments):
            assert np.allclose(seg[:, 1], [lo[j], hi[j]], rtol=0, atol=1e-14)
    fig.text(.35 / WIDTH, .65 / total_height,
             'Bars: accepted estimates; whiskers: 95% CIs (bulk: approximate, from authors’ p values).\n'
             '† Bulk Kcnc3 is the authors’ merged Kcnc3,Napsa locus. Verdict counts refer to negative estimates.',
             fontsize=8, color='#444444', va='center')
    assert len(artists) == 54
    for bar, row in artists:
        assert bar.get_height() == row['log2FC']
    return artists


def main():
    # Entry point: builds the standalone replication figure.
    import runpy
    runpy.run_path(str(Path(__file__).with_name('fig4_replication.py')), run_name='__main__')


if __name__ == '__main__':
    main()

#%% debug
