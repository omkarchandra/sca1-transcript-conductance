"""P2b (pre-specified before any number is seen): Purkinje fraction and the
predicted bulk fall..

f per gene = Purkinje counts / counts of all 13 cell types, from the S08 pseudobulk raw counts summed over the 12 WT mice of
the 24 discovery mice (all four ages); the nuclei per cell type behind those sums are reported. Predicted bulk log2FC =
f x our 12-week Purkinje log2FC (E3 TWELVE_WEEK_GENES; Kcnmb4 and Cacna2d2 are not in it and E3's LIMMA_RESULTS is not on
local disk, so they are left out). Observed = GSE122099 Atxn1 154Q/2Q KI cerebellum 12 weeks, Log2_Fold_Change as given
(Kcnc3 = the merged locus "Kcnc3,Napsa", flagged), cross-checked against P2's NINE_GENES.tsv. Spearman rho over the five
channel genes and over the available nine; rule: rho >= 0.8 over the five counts as support for the dilution explanation.
Read-only; writes a NEW real_v01/ next to this script, exclusive create.

  python -I -B 02_replication/p2_bulk/p2b_purkinje_fraction.py
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/p2b_purkinje_fraction'))
FIVE = ('Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g')
NINE = FIVE + ('Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2')
RULE_MIN = 0.8
S08 = _release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/')
P2 = _release_paths.expand('${DATA_ROOT}/publication/P2_external_replication/real_v01/')
GEO_TABLE = P2 + 'external_v02/GSE122099/files/GSE122099_Atxn1_154Q_KI_Cerebellum_12weeks.txt.gz'
PINS = {
    'raw_counts': (S08 + 'raw_counts.tsv.gz', '3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778'),
    'unit_metadata': (S08 + 'unit_metadata.tsv', '11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4'),
    'gene_metadata': (S08 + 'gene_metadata.tsv', '765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3'),
    'e3_twelve_week': (_release_paths.expand('${REPO_ROOT}/results/e3_s08_12w/TWELVE_WEEK_GENES.tsv'), '634d895c7a3df092e72ec0fe3eca8eab85c6dfd2ebcdbb33a03689e1cf56b9cc'),
    'geo_table': (GEO_TABLE, 'ad94cdcf8711aa3929d91d75baf7eff2a81f495a1aae1ff9c5f0756bd95d9c1e'),
    'p2_nine_genes': (P2 + 'tests/NINE_GENES.tsv', '22cba0dabd7e08b38b47cc4393ea77fbbc7aab96fd45fb375a00ffb8533cfd55'),
}
RULE_TEXT = ('P2b rule (pre-registered): rho >= 0.8 over the five counts as support for the dilution explanation; '
             'nothing else is inferred.')
NOTE = 'P2b: Purkinje fraction f (S08 pseudobulk, 12 WT discovery mice), predicted bulk log2FC = f x 12-week Purkinje log2FC (E3), observed GSE122099 154Q/2Q cerebellum 12 w.'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def wtsv(df, name):
    buf = io.StringIO()
    buf.write('# ' + NOTE + '\n')
    df.to_csv(buf, sep='\t', index=False, float_format='%.6g')
    with open(OUT / name, 'x') as f:
        f.write(buf.getvalue())


def rho(table, genes):
    t = table[table.gene.isin(genes)].dropna(subset=['predicted_bulk_log2fc', 'observed_log2fc'])
    if len(t) < 3:
        return dict(n=len(t), rho=float('nan'), genes=t.gene.tolist())
    return dict(n=len(t), rho=float(spearmanr(t.predicted_bulk_log2fc, t.observed_log2fc).statistic), genes=t.gene.tolist())


def main():
    P = {}
    for key, (rel, want) in PINS.items():
        if sha(ROOT / rel) != want:
            raise ValueError(f'input hash differs: {key} {rel}')
        P[key] = ROOT / rel
    raw = pd.read_csv(P['raw_counts'], sep='\t')
    genes = pd.read_csv(P['gene_metadata'], sep='\t')
    units = pd.read_csv(P['unit_metadata'], sep='\t', dtype={'animal_id': str})
    if not raw.stable_id.equals(genes.stable_id) or list(raw.columns[1:]) != units.unit_label.tolist():
        raise ValueError('raw counts, gene metadata and unit labels do not line up')
    wt = units[units.genotype == 'WT']
    if wt.animal_id.nunique() != 12 or wt.cell_type_raw.nunique() != 13 or len(wt) != 12 * 13:
        raise ValueError('expected 12 WT mice x 13 cell types')
    sym = genes.set_index('dataset_symbol').stable_id
    if any(sym.index.tolist().count(g) != 1 for g in NINE):
        raise ValueError('a nine-gene symbol is not unique in gene_metadata')
    rows = raw.set_index('stable_id').loc[sym[list(NINE)].to_numpy()]
    rows.index = list(NINE)
    pk_cols = wt.loc[wt.cell_type_raw == 'Purkinje cell', 'unit_label'].tolist()
    pk_counts, all_counts = rows[pk_cols].sum(axis=1), rows[wt.unit_label.tolist()].sum(axis=1)
    f = pk_counts / all_counts
    nuclei = (wt.groupby('cell_type_raw').agg(n_units=('unit_label', 'size'), n_nuclei=('n_nuclei', 'sum')).reset_index()
              .sort_values('n_nuclei', ascending=False))
    nuclei['share_of_nuclei'] = nuclei.n_nuclei / nuclei.n_nuclei.sum()

    e3 = pd.read_csv(P['e3_twelve_week'], sep='\t')
    e3 = e3[(e3.population == 'Purkinje cell') & (e3.age_weeks == 12)].set_index('dataset_symbol')
    if e3.index.duplicated().any():
        raise ValueError('E3 table has duplicate Purkinje 12-week rows')
    geo = pd.read_csv(P['geo_table'], sep='\t')
    gcol = geo.columns[1]
    lfc = pd.to_numeric(geo.Log2_Fold_Change, errors='coerce')
    p2 = pd.read_csv(P['p2_nine_genes'], sep='\t', comment='#')
    p2 = p2[(p2.comparison == 'GSE122099_Atxn1_154Q_KI_Cerebellum_12weeks') & (p2.kcnc3 == 'merged')].set_index('gene')
    table, dev = [], 0.0
    for g in NINE:
        hit = np.flatnonzero((geo[gcol] == g).to_numpy())
        flag = ''
        if len(hit) == 0 and g == 'Kcnc3':
            hit = np.flatnonzero(geo[gcol].astype(str).str.split(',').apply(lambda s: g in s).to_numpy())
            flag = 'observed from merged locus ' + str(geo[gcol].iloc[hit[0]]) if len(hit) == 1 else ''
        if len(hit) != 1:
            raise ValueError(f'observed row for {g} is not unique: {len(hit)}')
        obs = float(lfc.iloc[hit[0]])
        dev = max(dev, abs(obs - float(p2.loc[g, 'log2fc_mut_vs_wt'])))
        in_e3 = g in e3.index
        pk = float(e3.loc[g, 'log2FC']) if in_e3 else np.nan
        table.append(dict(gene=g, gene_class='channel' if g in FIVE else 'other', purkinje_counts_WT=int(pk_counts[g]), all_counts_WT=int(all_counts[g]),
                          f_purkinje=float(f[g]), purkinje_12w_log2fc=pk, purkinje_12w_CI_L=float(e3.loc[g, 'CI_L']) if in_e3 else np.nan,
                          purkinje_12w_CI_R=float(e3.loc[g, 'CI_R']) if in_e3 else np.nan, purkinje_12w_p=float(e3.loc[g, 'p']) if in_e3 else np.nan,
                          predicted_bulk_log2fc=float(f[g]) * pk if in_e3 else np.nan, observed_log2fc=obs,
                          observed_status=str(geo.Test_Status.iloc[hit[0]]), observed_p=float(pd.to_numeric(geo.p_value, errors='coerce').iloc[hit[0]]),
                          note='; '.join(x for x in (flag, '' if in_e3 else 'not in E3 TWELVE_WEEK_GENES; E3 LIMMA_RESULTS not on local disk; left out') if x)))
    if dev > 1e-9:
        raise ValueError(f'observed log2FC differs from P2 NINE_GENES.tsv by {dev}')
    table = pd.DataFrame(table)
    five, nine = rho(table, FIVE), rho(table, NINE)
    five_wo_kcnc3 = rho(table, tuple(g for g in FIVE if g != 'Kcnc3'))
    support = bool(np.isfinite(five['rho']) and five['n'] == 5 and five['rho'] >= RULE_MIN)

    lines = ['# P2b Purkinje fraction and the predicted bulk fall', '', NOTE, '',
             f'Mechanical reading: Spearman rho over the five channel genes = {five["rho"]:.3f} (n {five["n"]}); '
             f'rule rho >= {RULE_MIN}: **{"support" if support else "no support"}** for the dilution explanation.',
             f'Over the available nine: rho = {nine["rho"]:.3f} (n {nine["n"]}: {", ".join(nine["genes"])}). '
             f'Without Kcnc3 (merged locus): rho over four = {five_wo_kcnc3["rho"]:.3f}.', '',
             '| gene | f (Purkinje / all, WT) | Purkinje 12 w log2FC (E3) | predicted bulk | observed GSE122099 154Q/2Q CB 12 w | note |', '|---|---|---|---|---|---|']
    for r in table.itertuples():
        lines.append(f'| {r.gene} | {r.f_purkinje:.3f} | {r.purkinje_12w_log2fc:.3f} | {r.predicted_bulk_log2fc:.3f} | {r.observed_log2fc:.3f} | {r.note} |')
    lines += ['', 'Nuclei behind the pseudobulk sums (12 WT mice):', '', '| cell type | units | nuclei | share |', '|---|---|---|---|']
    lines += [f'| {r.cell_type_raw} | {r.n_units} | {r.n_nuclei} | {r.share_of_nuclei:.3f} |' for r in nuclei.itertuples()]
    lines += ['', RULE_TEXT, '']

    OUT.mkdir(parents=False, exist_ok=False)
    wtsv(table, 'P2b_TABLE.tsv')
    wtsv(nuclei, 'NUCLEI_PER_TYPE.tsv')
    with open(OUT / 'SUMMARY.md', 'x') as f_:
        f_.write('\n'.join(lines))
    rec = dict(note=NOTE, task='P2b', script_sha256=sha(Path(__file__)),
               inputs={rel: h for rel, h in PINS.values()}, rule=RULE_TEXT, rho_five=five, rho_nine_available=nine,
               rho_five_without_kcnc3=five_wo_kcnc3, support_mechanical=support, observed_vs_p2_max_abs=dev,
               f_definition='sum over the 12 WT discovery mice of Purkinje pseudobulk raw counts / sum over the same mice of all 13 cell types',
               predicted_definition='f x E3 12-week Purkinje log2FC (KI - WT)', left_out=['Kcnmb4', 'Cacna2d2'],
               left_out_reason='not in E3 TWELVE_WEEK_GENES; E3 LIMMA_RESULTS (eab5ccbe) is on the cluster only, not on local disk',
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
               finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f_:
        f_.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f_:
        for p in files:
            f_.write(f'{sha(p)}  {p.name}\n')
    print('\n'.join(lines))
    print(json.dumps({p.name: sha(p) for p in files}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
