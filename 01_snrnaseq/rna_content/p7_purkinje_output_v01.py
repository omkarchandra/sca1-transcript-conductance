"""P7: transcriptional output per Purkinje nucleus.

From the annotated h5ad layers/counts (all genes, all cell types): per nucleus total UMI (row sum) and detected genes (stored
counts > 0). Nuclei: the 280,295 included discovery nuclei (24 mice) and the 38,017 nuclei of the six 12-week mice (excluded
only for age), cell types from the packet00 nucleus manifest (S06 annotation). Genotypes: s06 unit_metadata for the 24 mice,
the packet00 animal_manifest for the six 12-week mice. Per mouse: medians for Purkinje, granule and all non-Purkinje nuclei;
R_gran and R_all for UMI and for genes; Purkinje share. Per age: KI/WT = mean KI per-mouse R / mean WT per-mouse R.
Rule (verbatim in RULE_TEXT) applied mechanically. Gate: Purkinje discovery totals equal S07 full_count_total exactly and
nuclei per unit equal unit_metadata. Read-only; writes a NEW real_v01/ next to this script, exclusive create.

  .venv-packet00/bin/python3 -I -B results/publication_20261026/P7_purkinje_output/p7_purkinje_output_v01.py
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import sys

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

import importlib.util as _release_importlib; _release_config = next(p / "config" / "paths.py" for p in Path(__file__).resolve().parents if (p / "config" / "paths.py").is_file()); _release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config); _release_paths = _release_importlib.module_from_spec(_release_spec); _release_spec.loader.exec_module(_release_paths); ROOT = Path(_release_paths.expand('${DATA_ROOT}'))
OUT = Path(__file__).resolve().parent / 'real_v01'
AGES = (5, 12, 18, 24, 30)
W12 = {'2061': 'WT', '2062': 'WT', '22018': 'WT', '2063': 'KI', '2065': 'KI', '22019': 'KI'}
PINS = {
    'h5ad': ('data/adata_leiden_0.175_annotated_cell_types_v3.h5ad', '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'),
    'unit_meta': ('results/s06_stageA_20260930_v01/unit_metadata.tsv', 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657'),
    'animals': ('results/packet00_sca1_audit_20260929_run01/animal_manifest.tsv', '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'),
    'nucleus_manifest': ('results/packet00_sca1_audit_20260929_run01/nucleus_manifest.tsv.gz', None),
    's07_purkinje_meta': ('results/s07_core_representation_20260930_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz', '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721'),
}
RULE_TEXT = ('P7 rule (verbatim): at 12 weeks or later, if KI/WT R_gran < 0.85 and all three KI mice are below the '
             'lowest WT mouse of that age, report “candidate reduction of transcriptional output per Purkinje nucleus”. Otherwise '
             'report “no evidence”. 5 weeks is descriptive. Caveats to state: nuclear UMI reflects nuclear transcript content and '
             'capture, not membrane area; the Purkinje share reflects dissociation as well as cell loss.')
NOTE = 'P7: per-mouse Purkinje / granule / all-other median UMI and detected genes from the h5ad raw counts; KI/WT per age.'
R_COLS = ('R_gran_umi', 'R_all_umi', 'R_gran_genes', 'R_all_genes', 'purkinje_share')


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


def row_stats(g, n_rows, chunk=4000):
    """Per-row sum of stored counts and number of stored counts > 0, read in row chunks."""
    ptr = g['indptr'][...].astype(np.int64)
    totals, detected = np.zeros(n_rows), np.zeros(n_rows, dtype=np.int64)
    for s in range(0, n_rows, chunk):
        e = min(s + chunk, n_rows)
        a, b = int(ptr[s]), int(ptr[e])
        v = g['data'][a:b].astype(np.float64)
        if not (np.isfinite(v).all() and (v >= 0).all() and (v == np.floor(v)).all()):
            raise ValueError('counts are not non-negative integers')
        row = np.repeat(np.arange(e - s), np.diff(ptr[s:e + 1]))
        totals[s:e] = np.bincount(row, weights=v, minlength=e - s)
        detected[s:e] = np.bincount(row, weights=(v > 0), minlength=e - s).astype(np.int64)
    return totals, detected


def main():
    P = {}
    for key, (rel, want) in PINS.items():
        h = sha(ROOT / rel)
        if want is not None and h != want:
            raise ValueError(f'input hash differs: {key} {rel}')
        P[key] = (ROOT / rel, h)
    man = pd.read_csv(P['nucleus_manifest'][0], sep='\t', dtype={'animal_id': str}, low_memory=False,
                      usecols=['obs_name', 'animal_id', 'cell_type_raw', 'included', 'exclusion_reason'])
    inc = man.included.astype(str).str.lower().eq('true')
    keep = inc | (man.animal_id.isin(W12) & man.exclusion_reason.eq('age_weeks in design.excluded_ages'))
    um = pd.read_csv(P['unit_meta'][0], sep='\t', dtype={'animal_id': str})
    animals = pd.read_csv(P['animals'][0], sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    mice = um.drop_duplicates('animal_id')[['animal_id', 'age_weeks', 'genotype', 'sex']].copy()
    if len(mice) != 24:
        raise ValueError('expected 24 discovery mice')
    for m, gt in W12.items():
        if animals.loc[m, 'genotype'] != gt or int(animals.loc[m, 'age_weeks']) != 12:
            raise ValueError(f'12-week genotype differs: {m}')
    mice = pd.concat([mice, pd.DataFrame(dict(animal_id=list(W12), age_weeks=12, genotype=list(W12.values()),
                                              sex=[animals.loc[m, 'sex'] for m in W12]))], ignore_index=True)

    with h5py.File(P['h5ad'][0], 'r') as f:
        idx = read_elem(f['obs']).index.to_numpy(dtype=str)
        g = f['layers/counts']
        n_rows = int(g.attrs['shape'][0])
        if len(idx) != n_rows:
            raise ValueError('obs index length differs from the counts shape')
        totals, detected = row_stats(g, n_rows)
    stats = pd.DataFrame(dict(obs_name=idx, total_umi=totals, detected_genes=detected))
    if not stats.obs_name.is_unique:
        raise ValueError('obs names not unique')
    nuc = man[keep].merge(stats, on='obs_name', how='left', validate='one_to_one')
    if nuc.total_umi.isna().any() or len(nuc) != 280295 + 38017:
        raise ValueError('manifest nuclei missing from the h5ad or wrong count')

    # gates: Purkinje discovery totals vs S07; nuclei per discovery unit vs unit_metadata
    s07 = pd.read_csv(P['s07_purkinje_meta'][0], sep='\t').merge(stats, on='obs_name', how='left')
    if not np.array_equal(s07.full_count_total.to_numpy(dtype=float), s07.total_umi.to_numpy()):
        raise ValueError('gate: Purkinje totals differ from S07 full_count_total')
    per_unit = nuc[inc[keep].to_numpy()].groupby(['animal_id', 'cell_type_raw']).size()
    want = um.set_index(['animal_id', 'cell_type_raw']).n_nuclei
    if not per_unit.reindex(want.index).fillna(0).astype(int).equals(want.astype(int)):
        raise ValueError('gate: nuclei per unit differ from unit_metadata')

    rows = []
    for m in mice.itertuples():
        x = nuc[nuc.animal_id == m.animal_id]
        pk, gr, ot = x[x.cell_type_raw == 'Purkinje cell'], x[x.cell_type_raw == 'Granule cell'], x[x.cell_type_raw != 'Purkinje cell']
        r = dict(animal_id=m.animal_id, age_weeks=int(m.age_weeks), genotype=m.genotype, sex=m.sex, n_nuclei=len(x), n_purkinje=len(pk),
                 n_granule=len(gr), purkinje_median_umi=float(pk.total_umi.median()), granule_median_umi=float(gr.total_umi.median()),
                 other_median_umi=float(ot.total_umi.median()), purkinje_median_genes=float(pk.detected_genes.median()),
                 granule_median_genes=float(gr.detected_genes.median()), other_median_genes=float(ot.detected_genes.median()))
        r.update(R_gran_umi=r['purkinje_median_umi'] / r['granule_median_umi'], R_all_umi=r['purkinje_median_umi'] / r['other_median_umi'],
                 R_gran_genes=r['purkinje_median_genes'] / r['granule_median_genes'], R_all_genes=r['purkinje_median_genes'] / r['other_median_genes'],
                 purkinje_share=len(pk) / len(x))
        rows.append(r)
    per_mouse = pd.DataFrame(rows).sort_values(['age_weeks', 'genotype', 'animal_id']).reset_index(drop=True)

    ages, lines = [], ['# P7 transcriptional output per Purkinje nucleus', '', NOTE, '',
                       '| age | KI/WT R_gran (UMI) | KI per-mouse R_gran | WT per-mouse R_gran | all KI below lowest WT | KI/WT R_all (UMI) | KI/WT R_gran (genes) | KI/WT R_all (genes) | KI/WT Purkinje share | reading |',
                       '|---|---|---|---|---|---|---|---|---|---|']
    for a in AGES:
        d = per_mouse[per_mouse.age_weeks == a]
        ki, wt = d[d.genotype == 'KI'], d[d.genotype == 'WT']
        if len(ki) != 3 or len(wt) != 3:
            raise ValueError(f'expected 3 KI and 3 WT mice at {a} w')
        r = dict(age_weeks=a, **{f'KI_over_WT_{c}': float(ki[c].mean() / wt[c].mean()) for c in R_COLS})
        r['all_KI_below_lowest_WT_R_gran'] = bool(ki.R_gran_umi.max() < wt.R_gran_umi.min())
        if a >= 12:
            r['reading'] = ('candidate reduction of transcriptional output per Purkinje nucleus'
                            if r['KI_over_WT_R_gran_umi'] < 0.85 and r['all_KI_below_lowest_WT_R_gran'] else 'no evidence')
        else:
            r['reading'] = 'descriptive (5 weeks)'
        ages.append(r)
        lines.append(f'| {a} w | {r["KI_over_WT_R_gran_umi"]:.3f} | {", ".join(f"{v:.3f}" for v in ki.R_gran_umi)} | {", ".join(f"{v:.3f}" for v in wt.R_gran_umi)} | '
                     f'{"yes" if r["all_KI_below_lowest_WT_R_gran"] else "no"} | {r["KI_over_WT_R_all_umi"]:.3f} | {r["KI_over_WT_R_gran_genes"]:.3f} | '
                     f'{r["KI_over_WT_R_all_genes"]:.3f} | {r["KI_over_WT_purkinje_share"]:.3f} | {r["reading"]} |')
    ages = pd.DataFrame(ages)
    lines += ['', 'Caveats: nuclear UMI reflects nuclear transcript content and capture, not membrane area; the Purkinje share '
                  'reflects dissociation as well as cell loss.', '', RULE_TEXT, '']

    OUT.mkdir(parents=False, exist_ok=False)
    wtsv(per_mouse, 'P7_PER_MOUSE.tsv')
    wtsv(ages, 'P7_BY_AGE.tsv')
    with open(OUT / 'SUMMARY.md', 'x') as f:
        f.write('\n'.join(lines))
    rec = dict(note=NOTE, task='P7', script_sha256=sha(Path(__file__)),
               inputs={str(p.relative_to(ROOT)): h for p, h in P.values()},
               rule=RULE_TEXT, n_nuclei=dict(discovery=int(inc.sum()), w12=int((keep & ~inc).sum())),
               definitions=dict(total_umi='row sum of layers/counts over all 26,374 genes', detected_genes='number of stored counts > 0 in the row',
                                other='all nuclei of the mouse that are not Purkinje (granule included)', median='pandas median per mouse',
                                KI_over_WT='mean of KI per-mouse values / mean of WT per-mouse values'),
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, h5py=h5py.__version__,
               finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1, ensure_ascii=False) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print('\n'.join(lines))
    print(per_mouse[['animal_id', 'age_weeks', 'genotype', 'n_purkinje', 'purkinje_median_umi', 'granule_median_umi', 'R_gran_umi']].to_string(index=False))
    print(json.dumps({p.name: sha(p) for p in files}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
