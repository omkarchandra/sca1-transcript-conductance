"""A2: P1 repeated with negative-binomial count splitting..

A copy of p1_countsplit.py (e1298818) with one change of method: the thinning is the countsplit NB scheme. Per gene, the
NB size b = 1 / edgeR tagwise dispersion estimated once on the S08 Purkinje pseudobulk (a2_dispersion.R, real_v01/
dispersion/); every stored count c is split train ~ Binomial(c, p) with p ~ Beta(eps b, (1 - eps) b) (the beta-binomial),
test = c - train, eps = 0.5; the 12-week sensitivity uses eps = 0.25 the same way. A gene without an estimate (no count in
the Purkinje pseudobulk; only possible in the 12-week rows) is thinned binomially and counted. Gate added: the per-mouse
pseudobulk of the count matrix equals the S08 raw counts the dispersions came from. Outputs go to real_v01/splits/, with
P1's Poisson values beside. P1's description follows.

P1: Purkinje lane occupancy under count splitting..

Count splitting (Neufeld et al. 2024, Poisson thinning, eps = 0.5). For k = 0..9, every raw UMI count c of every gene of the
1,971 discovery Purkinje nuclei is split train ~ Binomial(c, 0.5), test = c - train: numpy default_rng(20261026 + k), one
binomial call over the CSR data array (rows in S07 nucleus_metadata order, elements in stored h5ad order).
  train  S06 normalize (each half's own whole-transcriptome totals); S06 fit_transform on the per-mouse moments (equal-mouse
         weights) for mean and scale on the frozen 1,120 S07 eligible genes; S07 weighted_pca, 30 PCs. Lane rule on the
         per-mouse means: a candidate PC separates KI from WT completely with the same sign at >= 3 of 4 ages; an age's lane
         PCs are the candidates completely separated at that age. Direction per age = unit vector of mean(KI per-mouse
         means) - mean(WT per-mouse means) on that age's lane PCs.
  test   normalised the same way, projected with the train mean, scale and rotation (${REPO_ROOT}/01_snrnaseq/s07_pca/s07_prepare_core_pca.py
         line 55 formula), scored on the train direction; threshold = linear p90 of the pooled WT test nuclei of the age;
         share = fraction strictly beyond it (pooled over nuclei, and the mean of per-mouse shares).
  12 w   the six held-out mice, full counts, projected with each split's train model and scored on the 5-w and 18-w train
         directions against that split's test-half thresholds (optional thinned sensitivity: W12_THIN).
Gates, hard stops before any split: raw totals, S06 normalisation, S06 transform and S07 PCA reproduced on the full counts;
lane sets, directions, WT p90 and per-mouse counts equal far-set real_v02; the 12-w full-count projection equals phase4_12w.
Imported unchanged: s06_features.normalize and fit_transform, s07_prepare_core_pca.weighted_pca. Read-only; writes a NEW
real_v01/ next to this script, exclusive create.

  python -I -B 01_snrnaseq/count_splitting/a2_nb_countsplit.py
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
import gzip
import hashlib
import io
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
import h5py

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / _release_paths.expand('${DATA_ROOT}/helpers/project_p00_v02r1_local1'))]   # as python -I -B
from anndata.io import read_elem
from s06_features import normalize, fit_transform
from s07_prepare_core_pca import weighted_pca

OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/a2_nb_countsplit'))
A2D = _release_paths.expand('${DATA_ROOT}/publication/A2_nb_countsplit/real_v01/dispersion/')
P1R = _release_paths.expand('${DATA_ROOT}/publication/P1_countsplit/real_v01/')
AGES = (5, 18, 24, 30)
RULE_AGES = (18, 24, 30)
SPLITS = tuple(range(10))
SEED0, EPS = 20261026, 0.5
W12_THIN = 0.25             # review : add the labelled Binomial(c, 0.25) 12-w sensitivity (None = full counts only)
RULE_SHARE = 'pooled'       # which KI share the reading rule reads: 'pooled' or 'per_mouse_mean' (both are written)
NO_LANE = 'exclude'         # a split without a lane at an age: 'exclude' from the mean, or count as share 'zero'
P90, N_PCS, RULE_MIN = 0.90, 30, 0.50
W12_MICE = {'2061': 'WT', '2062': 'WT', '22018': 'WT', '2063': 'KI', '2065': 'KI', '22019': 'KI'}
EXPECTED_LANES = {5: (5, 6, 11), 18: (5, 6, 11, 30), 24: (5, 6, 11, 30), 30: (5, 6, 30)}
BLOCK_FULL = {5: 0.54, 18: 0.95, 24: 0.75, 30: 0.91}
TOL = dict(exact=1e-12, eig_rel=1e-8, cos=1e-8, pcs=1e-6, far=1e-6)
S06A = _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/')
S07P = _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/')
FAR = _release_paths.expand('${REPO_ROOT}/results/far_set_ratios/real_v02/')
W12 = _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/')
PINS = {
    'h5ad': (_release_paths.expand('${DATA_ROOT}/data/adata_leiden_0.175_annotated_cell_types_v3.h5ad'), '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'),
    'unit_meta': (S06A + 'unit_metadata.tsv', 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657'),
    'cell_meta': (S06A + 'cell_metadata.tsv.gz', 'ee6c4bca63f237825d7d80389e759263dbb45477e8535abb015b350f0595dbe6'),
    'normalized': (S06A + 'normalized_log10k_features.npz', '082f1639be27db8791529bcaae7710b9cca4f8bf85b3f563ddd31de9455d9507'),
    'transform_keys': (S06A + 'transform_keys.json', 'a9f42ea00582c9e9174c732cb3c9ee9090507bbcad53cd442d23c3d5f28889b8'),
    'transform_params': (S06A + 'transform_parameters.npz', '0e4b000d17a09e33b9277f8e19e00eabf9e8dab319cbd71e179ad4b42cdb3b3b'),
    'features': (S06A + 'feature_registry.tsv', 'd49b138d3a806db64fb37b6994a30034569cd530bfe3b6a82becf10f367f8c24'),
    'pca_model': (S07P + 'pca_model.npz', 'eaf1468f99ed7935acf94a1e58b1a45f1fe1333ae6cc541c6afc550e4d128284'),
    'nucleus_pcs': (S07P + 'nucleus_PCs.npy', 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b'),
    'nucleus_meta': (S07P + 'nucleus_metadata.tsv.gz', '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721'),
    'pc_check': (_release_paths.expand('${REPO_ROOT}/results/s12_lanes/pc_check_by_pc.tsv'), '27701bd432255403c50f399ad9a9f98a374b8cdccdfb2e8c22e87849b3bc4b31'),
    'far_record': (FAR + 'RECORD.json', '7ed5c0c184ca977afe71476e6e92afbcd2a75b7e3e1f94536687b64d13df0ea4'),
    'far_occupancy': (FAR + 'OCCUPANCY.tsv', '2823a594148d1b625a4938d3cd8ecb12112a86478fd0b0535b1cf2bcaebf9b1e'),
    'w12_cells': (W12 + 'cell_metadata_12w.tsv.gz', 'f2d3bf5a6987bcd27a2006791e1e2f02d63b6fc9cbb79979b3ac66c4807200fd'),
    'w12_normalized': (W12 + 'normalized_log10k_features_12w.npz', '2fcbe32f4e458c0ef7133fbd842a4afd7e1b0f7b4e3b3e9e1dd4c004fa32820b'),
    'w12_pcs': (W12 + 'Purkinje_cell_PCs_12w.npy', '3c2eb711b5f3c6777c2f7deb43789337633132194bb09d168544d312ddcfb310'),
    'animals': (_release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'), '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'),
    'code_s06_features': (_release_paths.expand('${DATA_ROOT}/helpers/s06_features.py'), '935dc821d48b12f8a4ee77d65ebe8cdd38fe8b6477c39de91aad0c84c25a5290'),
    'code_s07_core_pca': (_release_paths.expand('${REPO_ROOT}/01_snrnaseq/s07_pca/s07_prepare_core_pca.py'), '82cd4e996a5b98ed4396619c7cf11f62fd89fccf70044ac7b5d3b23db3de9221'),
    'code_s06_build': (_release_paths.expand('${DATA_ROOT}/helpers/s06_build.py'), '072df3e7e03a345fa69f04bcc087bc13850762677390ac644d92ae87643f64d4'),
    'dispersion': (A2D + 'DISPERSION.tsv', '67908c8a4213cc481f5d4be11a83eeb9de3999bdc630142182ad87854e7759eb'),
    's08_raw_counts': (_release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/raw_counts.tsv.gz'), '3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778'),
    's08_unit_meta': (_release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/unit_metadata.tsv'), '11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4'),
    'p1_summary': (P1R + 'P1_SUMMARY.tsv', 'f7bf2f9dfc18ab4d802f76a718425eeb277a728cb927269a00c11763409fb343'),
    'p1_w12_summary': (P1R + 'P1_W12_SUMMARY.tsv', 'c42728d582ce3ccaeb82becfcf2c8a8fd979b21017a76df048db2dc88f47afdc'),
    'code_p00_counts_scan': (_release_paths.expand('${DATA_ROOT}/helpers/project_p00_v02r1_local1/p00/counts_scan.py'), '6781f28fb5e529318cf2a53d6a31e24558ed31bd0024c0ad4807aff2869f9c93'),
}
RULE_TEXT = ('P1 reading rule (pre-registered): The occupancy claim survives if the mean KI share over '
             'the 10 splits is ≥ 0.50 at 18, 24 and 30 weeks (WT share is 0.10 by construction); 5 weeks is reported '
             'descriptively. The number of splits in which no PC passes the lane rule at an age is reported; more than 5 of 10 '
             'at any age means the lane rule is reported as unstable at that age. The full-data values (0.54 / 0.95 / 0.75 / '
             '0.91) are shown beside the split values; a drop is expected because half the counts are used. Caveat to state: '
             'Poisson thinning ignores overdispersion; a negative-binomial variant is a follow-up if a check asks.')
NOTE = 'A2 count splitting (negative binomial: beta-binomial split, edgeR pseudobulk dispersion, eps 0.5, 10 splits); lane PCs re-selected on each train half; shares on the test half.'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def sha_array(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def read_rows(g, ptr, positions, n_genes):
    """Raw counts of the given h5ad rows, every gene, rows in the given order, elements in stored order."""
    data, indices, indptr = [], [], [0]
    for p in positions:
        a, b = int(ptr[p]), int(ptr[p + 1])
        data.append(g['data'][a:b])
        indices.append(g['indices'][a:b])
        indptr.append(indptr[-1] + (b - a))
    value = np.concatenate(data).astype(np.float64)
    if not (np.isfinite(value).all() and (value >= 0).all() and (value == np.floor(value)).all()):
        raise ValueError('h5ad counts are not non-negative integers')
    return sp.csr_matrix((value.astype(np.int64), np.concatenate(indices).astype(np.int64), np.asarray(indptr, dtype=np.int64)),
                         shape=(len(positions), n_genes))


def thin(counts, rng, eps, b_col):
    """Beta-binomial split of every stored count: p ~ Beta(eps b, (1 - eps) b) for the elements whose gene has an NB size b
    (one beta call over them), p = eps otherwise, then one binomial call over the data array; test = c - train."""
    b = b_col[counts.indices]
    fin = np.isfinite(b)
    p = np.full(counts.data.shape, eps, dtype=np.float64)
    p[fin] = rng.beta(eps * b[fin], (1 - eps) * b[fin])
    train = rng.binomial(counts.data, p)
    test = counts.data - train
    if not ((train >= 0).all() and (test >= 0).all() and np.array_equal(train + test, counts.data)):
        raise ValueError('thinning: train + test differs from the original counts')
    make = lambda d: sp.csr_matrix((d, counts.indices.copy(), counts.indptr.copy()), shape=counts.shape)
    tr, te = make(train), make(test)
    if (tr + te - counts).count_nonzero() != 0:
        raise ValueError('thinning: train + test matrix differs from the original matrix')
    return tr, te


def normalized_features(counts, positions):
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    return normalize(counts[:, positions], totals), totals


def unit_moments(norm, unit_row, n_units):  # ${DATA_ROOT}/helpers/s06_build.py lines 82-84, per Purkinje unit
    indicator = sp.csr_matrix((np.ones(len(unit_row)), (unit_row, np.arange(len(unit_row)))), shape=(n_units, len(unit_row)))
    return (indicator @ norm).toarray(), (indicator @ norm.power(2)).toarray(), (indicator @ (norm > 0)).toarray()


def fit_pca(norm, mean, scale, weights, eligible):  # ${REPO_ROOT}/01_snrnaseq/s07_pca/s07_prepare_core_pca.py lines 46-50
    matrix = norm[:, eligible]
    empirical = np.asarray(matrix.T @ weights).ravel()
    if not np.allclose(empirical, mean, atol=1e-10):
        raise ValueError('weighted mean differs from the transform mean')
    return weighted_pca(matrix, mean, scale, weights)


def project(norm, mean, scale, rotation, eligible):  # ${REPO_ROOT}/01_snrnaseq/s07_pca/s07_prepare_core_pca.py line 55
    return np.asarray(norm[:, eligible].multiply(1 / scale) @ rotation) - (mean / scale) @ rotation


def mouse_means(pcs, animal_of_row, mice):
    d = pd.DataFrame(pcs, columns=range(1, pcs.shape[1] + 1))
    d['animal_id'] = animal_of_row
    return d.groupby('animal_id').mean().loc[mice.animal_id]


def separation(pm, mice, n_pcs):
    calls = {}
    for a in AGES:
        ki = pm.loc[mice[(mice.age_weeks == a) & (mice.genotype == 'KI')].animal_id]
        wt = pm.loc[mice[(mice.age_weeks == a) & (mice.genotype == 'WT')].animal_id]
        for k in range(1, n_pcs + 1):
            calls[(k, a)] = ('KI_higher' if ki[k].min() > wt[k].max() else 'KI_lower' if ki[k].max() < wt[k].min() else 'none')
    return calls


def lane_sets(calls, n_pcs):
    cand, major = [], {}
    for k in range(1, n_pcs + 1):
        hi = sum(calls[(k, a)] == 'KI_higher' for a in AGES)
        lo = sum(calls[(k, a)] == 'KI_lower' for a in AGES)
        if max(hi, lo) >= 3:
            cand.append(k)
            major[k] = 'KI_higher' if hi >= 3 else 'KI_lower'
    lanes = {a: tuple(k for k in cand if calls[(k, a)] != 'none') for a in AGES}
    opposite = [(k, a) for k in cand for a in AGES if calls[(k, a)] not in ('none', major[k])]
    return tuple(cand), lanes, opposite


def direction(pm, mice, a, lanes):
    ki = mice[(mice.age_weeks == a) & (mice.genotype == 'KI')].animal_id
    wt = mice[(mice.age_weeks == a) & (mice.genotype == 'WT')].animal_id
    d = np.array([pm.loc[ki, k].mean() - pm.loc[wt, k].mean() for k in lanes])
    return d / np.linalg.norm(d)


def shares(scores, thr, animal, genotype):
    """Pooled and per-mouse shares strictly beyond thr, by genotype."""
    beyond = scores > thr
    out = {}
    for g in ('KI', 'WT'):
        sel = genotype == g
        per = pd.Series(beyond[sel]).groupby(animal[sel]).agg(['sum', 'size'])
        out[g] = dict(pooled=float(beyond[sel].mean()), per_mouse_mean=float((per['sum'] / per['size']).mean()),
                      n_nuclei=int(sel.sum()), per_mouse=per)
    return out


def fmt_dir(lanes, u):
    return ';'.join(f'PC{k}:{v:+.6f}' for k, v in zip(lanes, u))


def wtsv(df, name):
    buf = io.StringIO()
    buf.write('# ' + NOTE + '\n')
    df.to_csv(buf, sep='\t', index=False, float_format='%.6g')
    with open(OUT / name, 'x') as f:
        f.write(buf.getvalue())


def main():
    started = datetime.now(timezone.utc)
    for mod, key in (('s06_features', 'code_s06_features'), ('s07_prepare_core_pca', 'code_s07_core_pca'), ('s06_build', 'code_s06_build'),
                     ('p00.counts_scan', 'code_p00_counts_scan')):
        if Path(sys.modules[mod].__file__).resolve() != (ROOT / PINS[key][0]).resolve():
            raise ValueError(f'imported {mod} is not the pinned file')
    P = {}
    for key, (rel, want) in PINS.items():
        if sha(ROOT / rel) != want:
            raise ValueError(f'input hash differs: {key} {rel}')
        P[key] = ROOT / rel
    rec = dict(note=NOTE, task='A2; design of P1', script_sha256=sha(Path(__file__)),
               inputs={rel: h for rel, h in PINS.values()}, started_utc=started.isoformat(timespec='seconds'),
               python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__,
               h5py=h5py.__version__, platform=platform.platform(), reading_rule=RULE_TEXT,
               choices=dict(rule_share=RULE_SHARE, no_lane=NO_LANE, w12_thin=W12_THIN, sd_ddof=1, eps=EPS, seeds=[SEED0 + k for k in SPLITS],
                            thinning='numpy default_rng(20261026 + k): p = beta(0.5 b, 0.5 b) over the stored elements whose gene has b (b = 1 / edgeR tagwise dispersion), '
                                     'p = 0.5 elsewhere; then binomial(data, p): one call each over the int64 CSR data array of the '
                                     '1,971 x 26,374 raw count matrix, rows in S07 nucleus_metadata order, elements in stored h5ad order; '
                                     'test = data - train' + ('; then the 12-w rows from the same generator' if W12_THIN else ''),
                            train_transform='S06 fit_transform on the train half per-mouse moments; mean and scale taken on the frozen 1,120 genes',
                            lane_rule='candidate: complete KI/WT separation of per-mouse means with the same sign at >= 3 of 4 ages; '
                                      "an age's lane PCs: candidates completely separated at that age (either sign)",
                            threshold='numpy quantile 0.9 (linear) of the pooled WT test nuclei of the age; beyond = strictly greater'))

    # discovery mice and nuclei
    um = pd.read_csv(P['unit_meta'], sep='\t', dtype={'animal_id': str})
    pk = um[um.cell_type_raw == 'Purkinje cell'].reset_index(drop=True)
    if not (len(pk) == 24 and set(pk.age_weeks) == set(AGES) and (pk.groupby(['age_weeks', 'genotype']).size() == 3).all() and (pk.n_nuclei >= 20).all()):
        raise ValueError('unexpected Purkinje unit table')
    nm = pd.read_csv(P['nucleus_meta'], sep='\t')
    unit_pos = pd.Series(np.arange(24), index=pk.unit_index.to_numpy())
    if not nm.unit_index.isin(unit_pos.index).all() or not nm.obs_name.is_unique:
        raise ValueError('nucleus_metadata rows outside the 24 Purkinje units')
    unit_row = unit_pos[nm.unit_index.to_numpy()].to_numpy()
    if not np.array_equal(np.bincount(unit_row, minlength=24), pk.n_nuclei.to_numpy()):
        raise ValueError('nuclei per unit differ from unit_metadata')
    mice = pk[['animal_id', 'age_weeks', 'genotype', 'sex']]
    animal = pk.animal_id.to_numpy()[unit_row]
    age = pk.age_weeks.to_numpy()[unit_row]
    geno = pk.genotype.to_numpy()[unit_row]
    weights = 1 / (pk.n_nuclei.to_numpy()[unit_row] * 24)

    w12c = pd.read_csv(P['w12_cells'], sep='\t', dtype={'animal_id': str})
    w12_sel = np.flatnonzero(w12c.cell_type_raw.eq('Purkinje cell').to_numpy())
    w12 = w12c.iloc[w12_sel].reset_index(drop=True)
    animals = pd.read_csv(P['animals'], sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    if set(w12.animal_id) != set(W12_MICE) or any(animals.loc[m, 'genotype'] != g or int(animals.loc[m, 'age_weeks']) != 12 for m, g in W12_MICE.items()):
        raise ValueError('12-week mice or genotypes differ')
    w12_geno = w12.animal_id.map(W12_MICE).to_numpy()
    w12_sex = {m: animals.loc[m, 'sex'] for m in W12_MICE}

    # raw counts from the h5ad
    features = pd.read_csv(P['features'], sep='\t')
    positions = features.var_position.to_numpy(dtype=int)
    with h5py.File(P['h5ad'], 'r') as f:
        obs_index = read_elem(f['obs']).index
        var = read_elem(f['var'])
        if not np.array_equal(var.iloc[positions].gene_ids.to_numpy(), features.stable_id.to_numpy()):
            raise ValueError('feature registry does not match the h5ad var')
        lookup = pd.Series(np.arange(len(obs_index)), index=obs_index)
        if not lookup.index.is_unique:
            raise ValueError('h5ad obs names are not unique')
        g = f['layers/counts']
        ptr = g['indptr'][...].astype(np.int64)
        n_genes = int(g.attrs['shape'][1])
        counts = read_rows(g, ptr, lookup[nm.obs_name].to_numpy(), n_genes)
        counts12 = read_rows(g, ptr, lookup[w12.obs_name].to_numpy(), n_genes)
    rec['n_nuclei'] = dict(discovery=int(counts.shape[0]), w12=int(counts12.shape[0]), genes=n_genes, nnz=int(counts.nnz), nnz_w12=int(counts12.nnz))
    rec['count_matrix_sha256'] = dict(data_int64_le=sha_array(counts.data.astype('<i8')), indices_int64_le=sha_array(counts.indices.astype('<i8')),
                                      indptr_int64_le=sha_array(counts.indptr.astype('<i8')))

    # A2: NB size per gene, in h5ad column order; gate: per-mouse pseudobulk of this matrix equals the S08 raw counts
    disp = pd.read_csv(P['dispersion'], sep='\t', comment='#').set_index('stable_id')
    gene_ids = var.gene_ids.to_numpy()
    b_col = disp.nb_size_b.reindex(gene_ids).to_numpy(dtype=np.float64)
    s08 = pd.read_csv(P['s08_raw_counts'], sep='\t').set_index('stable_id')
    s08u = pd.read_csv(P['s08_unit_meta'], sep='\t', dtype={'animal_id': str})
    s08u = s08u[s08u.cell_type_raw == 'Purkinje cell'].set_index('animal_id')
    ind = sp.csr_matrix((np.ones(len(unit_row)), (unit_row, np.arange(len(unit_row)))), shape=(24, len(unit_row)))
    pb = np.asarray((ind @ counts).todense())
    want = s08.reindex(gene_ids)[s08u.loc[pk.animal_id, 'unit_label'].tolist()].to_numpy(dtype=np.float64).T
    if not np.array_equal(pb, want):
        raise ValueError('gate: per-mouse pseudobulk differs from the S08 raw counts the dispersions came from')
    no_b_disc = int((~np.isfinite(b_col[counts.indices])).sum())
    if no_b_disc:
        raise ValueError(f'{no_b_disc} discovery counts have no NB size')
    rec['nb'] = dict(source='edgeR estimateDisp, robust, S08 Purkinje pseudobulk, ~0+genotype_age; b = 1/tagwise dispersion',
                     genes_with_b=int(np.isfinite(b_col).sum()), b_quantiles={q: float(np.nanquantile(b_col, q)) for q in (0, .1, .5, .9, 1)},
                     w12_counts_without_b=int((~np.isfinite(b_col[counts12.indices])).sum()), w12_count_elements=int(counts12.nnz),
                     pseudobulk_gate='exact')

    # gate A: totals, S06 normalisation, S06 transform
    norm_full, totals = normalized_features(counts, positions)
    if not np.array_equal(totals, nm.full_count_total.to_numpy(dtype=np.float64)):
        raise ValueError('gate A: raw row sums differ from S07 full_count_total')
    cm = pd.read_csv(P['cell_meta'], sep='\t')
    s06_row = pd.Series(np.arange(len(cm)), index=cm.obs_name)[nm.obs_name].to_numpy()
    s06_norm = sp.load_npz(P['normalized']).tocsr()[s06_row]
    dev_norm = float(np.abs(norm_full.toarray() - s06_norm.toarray()).max())
    keys = json.loads(P['transform_keys'].read_text())
    kk = next(i for i, key in enumerate(keys) if key['mode'] == 'discovery' and key['population'] == 'Purkinje cell')
    params = np.load(P['transform_params'])
    elig = params['eligible'][kk]
    if int(elig.sum()) != 1120:
        raise ValueError('frozen eligible list is not 1,120 genes')
    t_full = fit_transform(*unit_moments(norm_full, unit_row, 24), pk.n_nuclei.to_numpy(), np.ones(24, dtype=bool))
    dev_mean = float(np.abs(t_full['mean'] - params['mean'][kk]).max())
    dev_scale = float(np.abs(t_full['scale'] - params['scale'][kk]).max())
    if dev_norm > TOL['exact'] or dev_mean > TOL['exact'] or dev_scale > TOL['exact'] or not np.array_equal(t_full['eligible'], elig):
        raise ValueError(f'gate A: S06 normalisation/transform not reproduced (norm {dev_norm}, mean {dev_mean}, scale {dev_scale})')

    # gate B: S07 PCA on the full counts
    pca = np.load(P['pca_model'])
    mean_f, scale_f = params['mean'][kk, elig], params['scale'][kk, elig]
    if not (np.array_equal(pca['feature_indices'], np.flatnonzero(elig)) and np.allclose(pca['mean'], mean_f) and np.allclose(pca['scale'], scale_f)):
        raise ValueError('gate B: frozen PCA model does not match the frozen transform')
    values, rot = fit_pca(norm_full, mean_f, scale_f, weights, elig)          # as S07: frozen mean and scale
    pcs_full = project(norm_full, mean_f, scale_f, rot, elig)
    saved_pcs = np.load(P['nucleus_pcs'])
    dev_eig = float(np.max(np.abs(values - pca['eigenvalues']) / pca['eigenvalues'])) if len(values) == N_PCS else np.inf
    cos = np.sum(rot * pca['rotation'], axis=0) / (np.linalg.norm(rot, axis=0) * np.linalg.norm(pca['rotation'], axis=0)) if rot.shape == pca['rotation'].shape else np.zeros(1)
    dev_pcs = float(np.abs(pcs_full - saved_pcs).max()) if pcs_full.shape == saved_pcs.shape else np.inf
    if dev_eig > TOL['eig_rel'] or float(1 - cos.min()) > TOL['cos'] or dev_pcs > TOL['pcs']:
        raise ValueError(f'gate B: S07 PCA not reproduced (eig {dev_eig}, 1-cos {1 - cos.min()}, pcs {dev_pcs})')

    # gate C: lane rule, directions, thresholds and per-mouse counts on the full data equal far-set real_v02
    pm_full = mouse_means(pcs_full, animal, mice)
    calls_full = separation(pm_full, mice, N_PCS)
    cand_full, lanes_full, opp_full = lane_sets(calls_full, N_PCS)
    pcc = pd.read_csv(P['pc_check'], sep='\t')
    pcc = pcc[pcc.population == 'Purkinje_cell'].set_index('PC')
    pc_check_disagree = sum(calls_full[(k, a)] != pcc.loc[k, f'sep_{a}W'] for k in range(1, N_PCS + 1) for a in AGES)
    if cand_full != (5, 6, 11, 30) or lanes_full != EXPECTED_LANES:
        raise ValueError(f'gate C: lane sets differ: {cand_full} {lanes_full}')
    far = json.loads(P['far_record'].read_text())
    occ = pd.read_csv(P['far_occupancy'], sep='\t', comment='#', dtype={'mouse': str})
    full, dev_far = {}, 0.0
    for a in AGES:
        u = direction(pm_full, mice, a, lanes_full[a])
        rows = age == a
        s = pcs_full[rows][:, [k - 1 for k in lanes_full[a]]] @ u
        thr = float(np.quantile(s[geno[rows] == 'WT'], P90))
        old_u = np.array([far['lane_unit_vector'][str(a)][f'PC{k}'] for k in lanes_full[a]])
        dev_far = max(dev_far, float(np.abs(old_u - u).max()), abs(far['wt_quantile_thresholds'][str(a)]['0.9'] - thr))
        sh = shares(s, thr, animal[rows], geno[rows])
        for gt in ('KI', 'WT'):
            for m, r in sh[gt]['per_mouse'].iterrows():
                n_old = int(occ.loc[(occ.mouse == m) & (occ.age == a), 'n_beyond_WT_p90'].iloc[0])
                if int(r['sum']) != n_old:
                    raise ValueError(f'gate C: n beyond p90 differs from far-set real_v02: {a} w {m} {int(r["sum"])} vs {n_old}')
        full[a] = dict(u=u, thr=thr, KI_pooled=sh['KI']['pooled'], KI_per_mouse_mean=sh['KI']['per_mouse_mean'],
                       WT_pooled=sh['WT']['pooled'], WT_per_mouse_mean=sh['WT']['per_mouse_mean'])
    if dev_far > TOL['far']:
        raise ValueError(f'gate C: direction or p90 differs from far-set real_v02 by {dev_far}')
    gene_dir_full = {a: (rot[:, [k - 1 for k in lanes_full[a]]] @ full[a]['u']) / scale_f for a in AGES}

    # gate D: 12-week full counts reproduce phase4_12w
    norm12, totals12 = normalized_features(counts12, positions)
    if not np.array_equal(totals12, w12.full_count_total.to_numpy(dtype=np.float64)):
        raise ValueError('gate D: 12-week row sums differ from phase4_12w full_count_total')
    dev12_norm = float(np.abs(norm12.toarray() - sp.load_npz(P['w12_normalized']).tocsr()[w12_sel].toarray()).max())
    dev12_pcs = float(np.abs(project(norm12, mean_f, scale_f, pca['rotation'], elig) - np.load(P['w12_pcs'])).max())
    if dev12_norm > TOL['exact'] or dev12_pcs > TOL['pcs']:
        raise ValueError(f'gate D: 12-week normalisation/projection not reproduced (norm {dev12_norm}, pcs {dev12_pcs})')
    rec['gates'] = dict(A=dict(norm_max_abs=dev_norm, mean_max_abs=dev_mean, scale_max_abs=dev_scale, eligible_equal=True),
                        B=dict(eigenvalue_max_rel=dev_eig, min_cosine=float(cos.min()), pcs_max_abs=dev_pcs),
                        C=dict(lane_candidates=list(cand_full), lanes={str(a): list(v) for a, v in lanes_full.items()},
                               opposite_sign=opp_full, direction_p90_max_abs_vs_far=dev_far, per_mouse_counts_equal=True,
                               pc_check_label_disagreements_of_120=int(pc_check_disagree)),
                        D=dict(norm_max_abs=dev12_norm, pcs_max_abs=dev12_pcs), tolerances=TOL)

    # the splits
    split_rows, mouse_rows, w12_rows, split_rec = [], [], [], []
    tt = pd.DataFrame(dict(obs_name=nm.obs_name, animal_id=animal, full_total=totals.astype(np.int64)))
    for k in SPLITS:
        rng = np.random.default_rng(SEED0 + k)
        tr, te = thin(counts, rng, EPS, b_col)
        versions = {'full_counts': norm12}
        if W12_THIN:
            tr12, _ = thin(counts12, rng, W12_THIN, b_col)
            versions[f'thinned_{W12_THIN}'] = normalized_features(tr12, positions)[0]
        norm_tr, tot_tr = normalized_features(tr, positions)
        norm_te, tot_te = normalized_features(te, positions)
        if not np.array_equal(tot_tr + tot_te, totals):
            raise ValueError('thinning: half totals do not add up')
        tt[f'train_total_k{k}'] = tot_tr.astype(np.int64)
        t = fit_transform(*unit_moments(norm_tr, unit_row, 24), pk.n_nuclei.to_numpy(), np.ones(24, dtype=bool))
        mean, scale = t['mean'][elig], t['scale'][elig]
        vals, R = fit_pca(norm_tr, mean, scale, weights, elig)
        n_pcs = R.shape[1]
        pcs_tr = project(norm_tr, mean, scale, R, elig)
        pcs_te = project(norm_te, mean, scale, R, elig)
        pm = mouse_means(pcs_tr, animal, mice)
        calls = separation(pm, mice, n_pcs)
        cand, lanes, opp = lane_sets(calls, n_pcs)
        split_rec.append(dict(split=k, seed=SEED0 + k, train_data_sha256=sha_array(tr.data.astype('<i8')), n_pcs=n_pcs,
                              eigenvalues_top5=[float(v) for v in vals[:5]], train_eligible_of_1120=int(t['eligible'][elig].sum()),
                              lane_candidates=list(cand), opposite_sign=opp))
        u_of, thr_of = {}, {}
        for a in AGES:
            row = dict(split=k, seed=SEED0 + k, age=a, n_pcs=n_pcs, lane_candidates=','.join(map(str, cand)),
                       lane_pcs=','.join(map(str, lanes[a])), n_lane_pcs=len(lanes[a]),
                       opposite_sign_pcs=','.join(str(p) for p, aa in opp if aa == a))
            if not lanes[a]:
                split_rows.append(row)
                continue
            u = direction(pm, mice, a, lanes[a])
            rows = age == a
            s = pcs_te[rows][:, [p - 1 for p in lanes[a]]] @ u
            thr = float(np.quantile(s[geno[rows] == 'WT'], P90))
            sh = shares(s, thr, animal[rows], geno[rows])
            gdir = (R[:, [p - 1 for p in lanes[a]]] @ u) / scale
            row.update(direction=fmt_dir(lanes[a], u), cos_gene_direction_vs_full=float(gdir @ gene_dir_full[a] / (np.linalg.norm(gdir) * np.linalg.norm(gene_dir_full[a]))),
                       threshold_p90=thr, n_WT_nuclei=sh['WT']['n_nuclei'], n_KI_nuclei=sh['KI']['n_nuclei'],
                       KI_share_pooled=sh['KI']['pooled'], KI_share_per_mouse_mean=sh['KI']['per_mouse_mean'],
                       WT_share_pooled=sh['WT']['pooled'], WT_share_per_mouse_mean=sh['WT']['per_mouse_mean'])
            split_rows.append(row)
            for gt in ('KI', 'WT'):
                for m, r in sh[gt]['per_mouse'].iterrows():
                    mouse_rows.append(dict(split=k, age=a, mouse=m, genotype=gt, sex=mice.set_index('animal_id').loc[m, 'sex'],
                                           n_nuclei=int(r['size']), n_beyond=int(r['sum']), share=float(r['sum'] / r['size'])))
            u_of[a], thr_of[a] = u, thr
        for version, n12 in versions.items():
            p12 = project(n12, mean, scale, R, elig)
            for ref in (5, 18):
                if ref not in u_of:
                    continue
                s12 = p12[:, [p - 1 for p in lanes[ref]]] @ u_of[ref]
                for m in W12_MICE:
                    sel = (w12.animal_id == m).to_numpy()
                    w12_rows.append(dict(split=k, version=version, reference_age=ref, mouse=m, genotype=W12_MICE[m], sex=w12_sex[m],
                                         n_nuclei=int(sel.sum()), threshold_p90=thr_of[ref], n_beyond=int((s12[sel] > thr_of[ref]).sum()),
                                         share=float((s12[sel] > thr_of[ref]).mean())))
        print(f'split {k}: lanes {dict((a, lanes[a]) for a in AGES)}', flush=True)

    splits = pd.DataFrame(split_rows)
    per_mouse = pd.DataFrame(mouse_rows)
    w12_tab = pd.DataFrame(w12_rows)
    summary = []
    for a in AGES:
        d = splits[splits.age == a]
        has = d[d.n_lane_pcs > 0]
        r = dict(age=a, n_splits=len(d), n_splits_with_lane=len(has), n_splits_no_lane=int((d.n_lane_pcs == 0).sum()),
                 lane_rule_unstable=bool((d.n_lane_pcs == 0).sum() > 5))
        for col in ('KI_share_pooled', 'KI_share_per_mouse_mean', 'WT_share_pooled', 'WT_share_per_mouse_mean', 'cos_gene_direction_vs_full'):
            v = has[col].to_numpy(dtype=float) if len(has) else np.array([])
            if NO_LANE == 'zero' and col.startswith('KI_share'):
                v = np.concatenate([v, np.zeros(len(d) - len(has))])
            r[f'{col}_mean'] = float(v.mean()) if v.size else float('nan')
            r[f'{col}_sd'] = float(v.std(ddof=1)) if v.size > 1 else float('nan')
            r[f'{col}_min'] = float(v.min()) if v.size else float('nan')
            r[f'{col}_max'] = float(v.max()) if v.size else float('nan')
        r.update(full_KI_pooled=full[a]['KI_pooled'], full_KI_per_mouse_mean=full[a]['KI_per_mouse_mean'],
                 full_WT_pooled=full[a]['WT_pooled'], full_block_value=BLOCK_FULL[a])
        summary.append(r)
    summary = pd.DataFrame(summary)
    p1 = pd.read_csv(P['p1_summary'], sep='\t', comment='#').set_index('age')
    summary['P1_KI_share_pooled_mean'] = summary.age.map(p1.KI_share_pooled_mean)
    summary['P1_KI_share_pooled_sd'] = summary.age.map(p1.KI_share_pooled_sd)
    summary['P1_n_splits_no_lane'] = summary.age.map(p1.n_splits_no_lane)
    rule_col = 'KI_share_pooled_mean' if RULE_SHARE == 'pooled' else 'KI_share_per_mouse_mean_mean'
    met = {a: bool(summary.set_index('age').loc[a, rule_col] >= RULE_MIN) for a in RULE_AGES}
    unstable = [a for a in AGES if summary.set_index('age').loc[a, 'lane_rule_unstable']]
    verdict = 'survives' if all(met.values()) else 'does not survive'
    if w12_tab.empty:
        raise ValueError('no 12-week scores: no split has a lane at 5 or 18 w')
    by = w12_tab.groupby(['version', 'reference_age', 'genotype', 'split'])
    per_split = pd.DataFrame(dict(pooled=by.n_beyond.sum() / by.n_nuclei.sum(), per_mouse_mean=by.share.mean()))
    w12_sum = per_split.groupby(level=[0, 1, 2]).agg(['mean', 'std', 'min', 'max'])    # std: ddof 1
    w12_sum.columns = [f'{a}_{b}' for a, b in w12_sum.columns]
    w12_sum = w12_sum.reset_index()

    lines = ['# A2 negative-binomial count-split occupancy (Purkinje lanes), beside P1 (Poisson)', '', NOTE, '',
             f'Mechanical reading ({RULE_SHARE} KI share, no-lane splits {NO_LANE}d): occupancy claim **{verdict}** '
             f'(mean KI share >= {RULE_MIN} at 18/24/30 w: ' + ', '.join(f'{a} w {"yes" if met[a] else "no"}' for a in RULE_AGES) + '); '
             f'lane rule unstable at: {", ".join(f"{a} w" for a in unstable) if unstable else "no age"}.', '',
             '| age | splits with lane | no-lane splits | KI share pooled, mean ± SD (min–max) | KI share per-mouse mean, mean ± SD | WT share pooled, mean | P1 Poisson KI pooled, mean ± SD | full data KI pooled / per-mouse mean / block |',
             '|---|---|---|---|---|---|---|---|']
    for r in summary.itertuples():
        lines.append(f'| {r.age} w | {r.n_splits_with_lane} | {r.n_splits_no_lane} | {r.KI_share_pooled_mean:.3f} ± {r.KI_share_pooled_sd:.3f} '
                     f'({r.KI_share_pooled_min:.3f}–{r.KI_share_pooled_max:.3f}) | {r.KI_share_per_mouse_mean_mean:.3f} ± {r.KI_share_per_mouse_mean_sd:.3f} | '
                     f'{r.WT_share_pooled_mean:.3f} | {r.P1_KI_share_pooled_mean:.3f} ± {r.P1_KI_share_pooled_sd:.3f} | {r.full_KI_pooled:.3f} / {r.full_KI_per_mouse_mean:.3f} / {r.full_block_value:.2f} |')
    lines += ['', '12-week held-out mice (share beyond the split\'s test-half WT p90; mean over splits):', '',
              '| version | reference age | genotype | pooled mean ± SD | per-mouse mean ± SD |', '|---|---|---|---|---|']
    for r in w12_sum.itertuples():
        lines.append(f'| {r.version} | {r.reference_age} w | {r.genotype} | {r.pooled_mean:.3f} ± {r.pooled_std:.3f} | {r.per_mouse_mean_mean:.3f} ± {r.per_mouse_mean_std:.3f} |')
    p1w = pd.read_csv(P['p1_w12_summary'], sep='\t', comment='#')
    lines += ['', 'P1 (Poisson) 12-week values for comparison: ' + '; '.join(f'{r.version} {r.reference_age} w {r.genotype} {r.pooled_mean:.3f}' for r in p1w.itertuples())]
    lines += ['', 'Method note: the NB size per gene comes from the between-mouse dispersion of the Purkinje pseudobulk (edgeR estimateDisp), '
              'the estimator the block allows; nucleus-level overdispersion (sctransform theta) is not installed here.',
              'Lane PC numbers are those of each split\'s own train PCA and are not the full-data PC numbers.', '', RULE_TEXT, '']

    OUT.mkdir(parents=False, exist_ok=False)
    wtsv(splits, 'A2_SPLITS.tsv')
    wtsv(per_mouse, 'A2_PER_MOUSE.tsv')
    wtsv(w12_tab, 'A2_W12.tsv')
    wtsv(summary, 'A2_SUMMARY.tsv')
    wtsv(w12_sum, 'A2_W12_SUMMARY.tsv')
    with open(OUT / 'A2_TRAIN_TOTALS.tsv.gz', 'xb') as f:
        f.write(gzip.compress(tt.to_csv(sep='\t', index=False).encode(), mtime=0))
    with open(OUT / 'SUMMARY.md', 'x') as f:
        f.write('\n'.join(lines))
    rec.update(splits=split_rec, verdict_mechanical=verdict, rule_ages_met={str(a): v for a, v in met.items()},
               lane_rule_unstable_ages=unstable, full_data={str(a): dict(direction=fmt_dir(lanes_full[a], full[a]['u']), p90=full[a]['thr'],
               KI_pooled=full[a]['KI_pooled'], KI_per_mouse_mean=full[a]['KI_per_mouse_mean'], WT_pooled=full[a]['WT_pooled'],
               WT_per_mouse_mean=full[a]['WT_per_mouse_mean']) for a in AGES},
               finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print('\n'.join(lines[:14]))
    print(json.dumps({p.name: sha(p) for p in files}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
