"""P2 step (a): the Purkinje control set, built and hashed
before any external file is opened..

The 200 genes most enriched in Purkinje nuclei in the discovery data. S08 pseudobulk raw counts: 26,374 genes x 312 units
(13 S06 cell types x 24 mice, WT and KI, all ages). Each unit's CPM is taken over all genes; score = log2((mean CPM of the 24
Purkinje units + 1) / (mean CPM of the 288 units of the other 12 cell types + 1)). Universe: the genes that pass the S08
expression filter in Purkinje (frozen registry). The nine genes are removed; top 200 by score, ties by stable_id.
Read-only; writes a NEW real_v01/control/ next to this script, exclusive create.

  python -I -B 02_replication/p2_bulk/p2_control_set.py
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

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/p2_control_set'))
N_CONTROL, PSEUDO_CPM = 200, 1.0
NINE = {'Kcnma1': 'ENSMUSG00000063142', 'Kcnc3': 'ENSMUSG00000062785', 'Kcnd3': 'ENSMUSG00000040896',
        'Cacna1a': 'ENSMUSG00000034656', 'Cacna1g': 'ENSMUSG00000020866', 'Itpr1': 'ENSMUSG00000030102',
        'Atp2a3': 'ENSMUSG00000020788', 'Kcnmb4': 'ENSMUSG00000054934', 'Cacna2d2': 'ENSMUSG00000010066'}
S08 = _release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/')
PINS = {
    'raw_counts': (S08 + 'raw_counts.tsv.gz', '3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778'),
    'gene_metadata': (S08 + 'gene_metadata.tsv', '765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3'),
    'unit_metadata': (S08 + 'unit_metadata.tsv', '11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4'),
    'registry': (_release_paths.expand('${DATA_ROOT}/results/s08_gene_registry_v01/PLANNED_GENE_TESTS.tsv.gz'), 'f648f3a4aca1a23cabde911b264a84401bd58588350ba6d558ace58142a53fd5'),
}
NOTE = ('P2 control set: top 200 Purkinje-enriched genes (S08 pseudobulk, all 24 discovery mice; log2((Purkinje mean CPM + 1) / '
        '(other 12 cell types mean CPM + 1)); S08 Purkinje expression-filter universe; nine genes removed). Built before any external file was opened.')


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
    with open(OUT / 'control' / name, 'x') as f:
        f.write(buf.getvalue())


def main():
    if OUT.exists():
        raise ValueError('real_v01 already exists')
    P = {}
    for key, (rel, want) in PINS.items():
        if sha(ROOT / rel) != want:
            raise ValueError(f'input hash differs: {key} {rel}')
        P[key] = ROOT / rel
    raw = pd.read_csv(P['raw_counts'], sep='\t')
    genes = pd.read_csv(P['gene_metadata'], sep='\t')
    units = pd.read_csv(P['unit_metadata'], sep='\t', dtype={'animal_id': str})
    if not raw.stable_id.equals(genes.stable_id) or len(raw) != 26374 or not genes.stable_id.is_unique:
        raise ValueError('raw counts and gene metadata differ')
    if list(raw.columns[1:]) != units.unit_label.tolist() or len(units) != 312:
        raise ValueError('raw count columns differ from the unit labels')
    per_type = units.groupby('cell_type_raw').size()
    if len(per_type) != 13 or not (per_type == 24).all():
        raise ValueError('expected 13 cell types x 24 units')
    for sym, sid in NINE.items():
        if genes.loc[genes.stable_id == sid, 'dataset_symbol'].tolist() != [sym]:
            raise ValueError(f'nine-gene identifier differs: {sym} {sid}')
    counts = raw.iloc[:, 1:].to_numpy(dtype=np.float64)
    lib = counts.sum(axis=0)
    cpm = counts / lib * 1e6
    pk = (units.cell_type_raw == 'Purkinje cell').to_numpy()
    pk_mean, other_mean = cpm[:, pk].mean(axis=1), cpm[:, ~pk].mean(axis=1)
    score = np.log2((pk_mean + PSEUDO_CPM) / (other_mean + PSEUDO_CPM))
    reg = pd.read_csv(P['registry'], sep='\t', usecols=['stable_id', 'population'])
    universe = set(reg.loc[reg.population == 'Purkinje cell', 'stable_id'])
    if not set(NINE.values()) <= universe:
        raise ValueError('a nine-gene id is outside the S08 Purkinje universe')
    tab = pd.DataFrame(dict(stable_id=genes.stable_id, symbol=genes.dataset_symbol, purkinje_mean_cpm=pk_mean,
                            other_mean_cpm=other_mean, log2_enrichment=score))
    tab = tab[tab.stable_id.isin(universe)].sort_values(['log2_enrichment', 'stable_id'], ascending=[False, True]).reset_index(drop=True)
    tab.insert(0, 'rank_in_universe', np.arange(1, len(tab) + 1))
    nine = tab[tab.stable_id.isin(NINE.values())].copy()
    control = tab[~tab.stable_id.isin(NINE.values())].head(N_CONTROL).copy()
    control.insert(0, 'control_rank', np.arange(1, N_CONTROL + 1))
    if len(control) != N_CONTROL or control.stable_id.duplicated().any():
        raise ValueError('control set is not 200 distinct genes')

    OUT.mkdir(parents=False, exist_ok=False)
    (OUT / 'control').mkdir()
    wtsv(control, 'CONTROL_SET.tsv')
    wtsv(nine, 'NINE_GENE_RANKS.tsv')
    rec = dict(note=NOTE, task='P2 step (a)', script_sha256=sha(Path(__file__)),
               inputs={rel: h for rel, h in PINS.values()}, built_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
               external_files_opened_before=False, n_universe=len(tab), n_control=N_CONTROL, pseudo_cpm=PSEUDO_CPM,
               excluded_nine=NINE, cutoff_log2_enrichment=float(control.log2_enrichment.iloc[-1]),
               control_min_purkinje_cpm=float(control.purkinje_mean_cpm.min()), python=sys.version.split()[0],
               numpy=np.__version__, pandas=pd.__version__,
               definition='per-unit CPM over all 26,374 genes; mean over the 24 Purkinje units vs mean over the 288 units of the '
                          'other 12 S06 cell types; all 24 discovery mice; score log2((Purkinje + 1) / (other + 1)); universe = '
                          'S08 Purkinje expression-filter genes; top 200 after removing the nine; ties by stable_id')
    with open(OUT / 'control' / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in (OUT / 'control').iterdir() if p.is_file())
    with open(OUT / 'control' / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print(control.head(15).to_string(index=False))
    print(nine.to_string(index=False))
    print(json.dumps({'n_universe': len(tab), 'cutoff': rec['cutoff_log2_enrichment'], 'min_pk_cpm': rec['control_min_purkinje_cpm'],
                      **{p.name: sha(p) for p in files}, 'SHA256SUMS': sha(OUT / 'control' / 'SHA256SUMS')}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
