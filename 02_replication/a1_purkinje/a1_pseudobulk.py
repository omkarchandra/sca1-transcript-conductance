"""v03: adds --root (the A1 package folder; default: this script's folder); no other change.
v02: reads only the Purkinje columns of each Borgenheimer table (v01, 51eefdad4c5d..., was OOM-killed on the laptop at 8.3 GiB).
A1 step (b): Purkinje pseudobulk per mouse or donor..

(i) GSE215336 (Borgenheimer 2022): per-sample cleaned_count_table (integer genes x nuclei) summed over the nuclei whose SingleR
label is "Neuron.Gad1Gad2.Pcp2"; genotype from the series matrix titles. (ii) GSE246183 human arm: raw/X integer counts summed
per donor over nuclei with sub4 == "PC" (libraries 409_2, 1516_2, 12_2 and 13_2 are summed with their donor); donors with fewer
than MIN_PC_HUMAN Purkinje nuclei are listed and left out (primary); all 20 donors as a labelled sensitivity. Writes counts (genes x samples) and sample tables into a NEW
real_v01/pseudobulk/; the limma step is a1_limma.R. Downloaded files are only read as tables and HDF5.

  python -I -B 02_replication/a1_purkinje/a1_pseudobulk.py [--root DIR]
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
import gzip
import hashlib
import json
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem

HERE = Path(next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == '--root'), Path(__file__).resolve().parent)).resolve()
EXT = HERE / 'real_v01' / 'external'
HUMAN = HERE / 'real_v01' / 'derived' / 'GSE246183_230502_human_ctype.h5ad'
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/a1_pseudobulk'))
EXT_SUMS = '93a23f1c4c6c9e4eefd5598e028974cc71805ef4af0ce7638b8908a14f0e5100'
MIN_PC_HUMAN = 10
B_PC_LABEL = 'Neuron.Gad1Gad2.Pcp2'
B_SAMPLES = {'GSM6632781_P7644': 'KI', 'GSM6632782_P7830': 'KI', 'GSM6632783_P7831': 'KI',
             'GSM6632784_P7643': 'WT', 'GSM6632785_P7832': 'WT', 'GSM6632786_P7833': 'WT'}   # series titles: SCA1 mouse 1-3, Wild-type 1-3
H_DONOR = {'14': 'Control_1', '409': 'Control_2', '4092': 'Control_2', '1516': 'Control_3', '15162': 'Control_3', '729': 'Control_4',
           '1432': 'Control_5', '1532': 'Control_6', 'T99': 'Control_7', 'T136': 'Control_8', 'T551': 'Control_9', 'T1952': 'Control_10',
           '11': 'SCA1_1', '12': 'SCA1_2', '122': 'SCA1_2', '13': 'SCA1_3', '132': 'SCA1_3', 'CB': 'SCA1_4', 'RT': 'SCA1_5', 'CT': 'SCA1_6',
           'UMB1354': 'SCA1_7', 'CNGA-0064': 'SCA1_8', 'CNGA-0047': 'SCA1_9', 'CNGA-0061': 'SCA1_10'}   # series-matrix descriptions (409_2 = 4092, ...)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def write(df, path):
    with open(path, 'x') as f:
        df.to_csv(f, sep='\t')


def main():
    if sha(EXT / 'SHA256SUMS') != EXT_SUMS:
        raise ValueError('external/SHA256SUMS differs')
    listed = {n: h for h, n in (l.split('  ', 1) for l in (EXT / 'SHA256SUMS').read_text().splitlines())}
    if OUT.exists():
        raise ValueError('pseudobulk/ already exists')
    rec = dict(task='A1 step (b), pseudobulk', script_sha256=sha(Path(__file__)), external_sums=EXT_SUMS, min_pc_human=MIN_PC_HUMAN)

    # (i) Borgenheimer
    raw = EXT / 'GSE215336' / 'unpacked_GSE215336_RAW'
    cols, meta, genes = {}, [], None
    for s, gt in B_SAMPLES.items():
        for part in ('celltyping', 'cleaned_count_table'):
            rel = f'GSE215336/unpacked_GSE215336_RAW/{s}_{part}.csv.gz'
            if sha(EXT / rel) != listed[rel]:
                raise ValueError(f'file differs from the manifest: {rel}')
        ct = pd.read_csv(raw / f'{s}_celltyping.csv.gz')
        pcb = ct.loc[ct.SingleR_Celltype == B_PC_LABEL, 'CellBarcode']
        header = pd.read_csv(raw / f'{s}_cleaned_count_table.csv.gz', nrows=0).columns
        if not pcb.isin(header).all():
            raise ValueError(f'{s}: Purkinje barcodes missing from the count table')
        # v02: read only the gene column and the Purkinje columns (v01 read the whole table densely and was OOM-killed)
        tab = pd.read_csv(raw / f'{s}_cleaned_count_table.csv.gz', usecols=[header[0]] + pcb.tolist(), index_col=0)
        v = tab.to_numpy()
        if not (np.issubdtype(v.dtype, np.integer) and (v >= 0).all()):
            raise ValueError(f'{s}: counts are not non-negative integers')
        genes = tab.index if genes is None else genes
        if not tab.index.equals(genes):
            raise ValueError(f'{s}: gene rows differ between samples')
        cols[s] = tab[pcb].sum(axis=1)
        meta.append(dict(sample=s, genotype=gt, age_weeks=12, n_nuclei=len(header) - 1, n_purkinje=len(pcb)))
    OUT.mkdir()
    (OUT / 'GSE215336').mkdir()
    write(pd.DataFrame(cols), OUT / 'GSE215336' / 'COUNTS.tsv')
    write(pd.DataFrame(meta).set_index('sample'), OUT / 'GSE215336' / 'SAMPLES.tsv')

    # (ii) human
    with h5py.File(HUMAN, 'r') as f:
        obs = read_elem(f['obs'])
        var = read_elem(f['raw/var'])
        g = f['raw/X']
        shape = tuple(int(x) for x in g.attrs['shape'])
        sel = np.flatnonzero((obs.sub4 == 'PC').to_numpy())
        ptr = g['indptr'][...].astype(np.int64)
        data, ind, iptr = [], [], [0]
        for r in sel:
            a, b = int(ptr[r]), int(ptr[r + 1])
            data.append(g['data'][a:b]); ind.append(g['indices'][a:b]); iptr.append(iptr[-1] + b - a)
    x = sp.csr_matrix((np.concatenate(data).astype(np.float64), np.concatenate(ind), np.array(iptr)), shape=(len(sel), shape[1]))
    if not (x.data == np.floor(x.data)).all():
        raise ValueError('human raw/X is not integer')
    pc = obs.iloc[sel]
    if set(pc.batch.astype(str)) - set(H_DONOR):
        raise ValueError('unmapped human batch')
    donor = pc.batch.astype(str).map(H_DONOR).to_numpy()
    donors = sorted(set(H_DONOR.values()))
    ind_m = sp.csr_matrix((np.ones(len(sel)), ([donors.index(d) for d in donor], np.arange(len(sel)))), shape=(len(donors), len(sel)))
    pb = pd.DataFrame(np.asarray((ind_m @ x).todense()).T, index=var.index.astype(str), columns=donors)
    gtypes = {d: ('KI' if d.startswith('SCA1') else 'WT') for d in donors}
    for d in donors:
        if set(obs.genotype[obs.batch.astype(str).map(H_DONOR) == d]) != {'SCA1' if gtypes[d] == 'KI' else 'Ctrl'}:
            raise ValueError(f'donor genotype differs from obs: {d}')
    n_pc = pd.Series(donor).value_counts().reindex(donors).fillna(0).astype(int)
    n_all = obs.batch.astype(str).map(H_DONOR).value_counts().reindex(donors).astype(int)
    hmeta = pd.DataFrame(dict(genotype=[gtypes[d] for d in donors], n_nuclei=n_all.to_numpy(), n_purkinje=n_pc.to_numpy(),
                              kept=(n_pc >= MIN_PC_HUMAN).to_numpy()), index=pd.Index(donors, name='sample'))
    hmeta['age_weeks'] = 0
    (OUT / 'GSE246183_human').mkdir()
    write(pb[hmeta.index[hmeta.kept]], OUT / 'GSE246183_human' / 'COUNTS.tsv')
    write(hmeta, OUT / 'GSE246183_human' / 'SAMPLES.tsv')
    (OUT / 'GSE246183_human_all20').mkdir()                    # review : labelled sensitivity, all 20 donors
    write(pb, OUT / 'GSE246183_human_all20' / 'COUNTS.tsv')
    write(hmeta.assign(kept=True), OUT / 'GSE246183_human_all20' / 'SAMPLES.tsv')
    rec['human'] = dict(n_purkinje_total=int(len(sel)), donors_kept=hmeta.index[hmeta.kept].tolist(), donors_left_out=hmeta.index[~hmeta.kept].tolist())
    rec['borgenheimer'] = dict(n_purkinje=[m['n_purkinje'] for m in meta], genes=int(len(genes)))
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in OUT.rglob('*') if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n')
    print(pd.DataFrame(meta).to_string(index=False)); print(hmeta.to_string()); print('SHA256SUMS', sha(OUT / 'SHA256SUMS'))


if __name__ == '__main__':
    main()

#%% debug
