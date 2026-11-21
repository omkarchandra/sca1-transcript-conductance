
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# S12 T6 , r2 (r1 used log-scale WT means; review's  block: markers on LINEAR-scale WT means, CP10k = expm1 of the S06 values): DIRECT CONTAMINATION SCORE of the scored populations Purkinje cell and Bergmann glia (pre-specified). DESCRIPTIVE, READ-ONLY, no interpretation, no refit.
# It is built on the T5 proofs (same inputs, same row mapping, same pins), reduced to what T6 needs.
# SCALE: the S06 matrix holds log1p(10000 * count / full_total) (s06_features.py normalize). Marker RATIOS use LINEAR-scale WT means: CP10k = expm1(value), averaged over WT nuclei; the log-scale pooled mean and ratio are written beside and are not used. The SCORE is the sum of the S06 log1p(CP10k) values (pre-specified).
# MARKERS (of each S07 population Q, from WT nuclei of the 24 included mice; nuclei = the S07 rule: units with n_nuclei >= 20): over all features of the S06 matrix, own = WT mean of linear CP10k expression in Q, max_other = max over the OTHER 12 populations of their WT means, ratio = own / max_other
#   (max_other = 0 and own > 0: ratio = +inf; own = 0: never a marker). Gene g is a marker of Q at threshold t when own > 0 and ratio >= t. Thresholds 2, 3, 4 (PRIMARY: review's proposal of 4x) and 6 are all scored in the one run; every gene's ratio is written, so any other reading of the ratio needs no rerun for the lists
#   (the scores for another threshold would). POOLED WT mean (all WT nuclei of the population) DEFINES THE MARKERS; the mouse-balanced WT mean (mean over WT mice of the per-mouse mean) and its ratio are written beside it. Genes that are markers of the scored population itself (same threshold) are removed from every Q's list, with the count of removed genes written.
# SCORE: per nucleus of the scored population P, the SUM of the normalized expression (S06 log1p CP10k values; no scaling, no centering) of the used markers of Q; per mouse the MEAN over the mouse's P nuclei; KI vs WT within age with the T3 exact test (exact_perm / separation copied from s12_t5_decomp_r1.py), per (P, Q, threshold).
#   A Q with no used marker at a threshold gets no score rows (status no_markers in marker_counts.tsv).
# PROOFS (STOP on failure; same as T5, quoted there): S07 RECORD -> S06 checksum chain; unit_index = arange; bincount(unit_index) = n_nuclei; matrix shape = (len(cell_metadata), len(registry)); units = the 24 included mice; cells.iloc[flatnonzero(p_mask[unit])] = the pinned nucleus_metadata of ALL 13 populations; and for the two scored populations the PC scores recomputed
#   in the S07 operation order from the matrix rows equal nucleus_PCs for all 30 PCs (|a-b| <= 1e-8 * max(1, |b|)), which proves that the matrix rows are aligned with cell_metadata.
# Writes: marker_genes.tsv (every gene x population: WT means, ratios, marker flags), marker_counts.tsv, marker_lists.tsv (the used markers), wt_denominators.tsv, contamination_mouse_means.tsv, contamination_kivswt.tsv, contamination_summary.tsv, checks.tsv, RESULT.json (nothing is overwritten; the output directory must be empty).
# usage: python -I -B s12_t6_contam_r2.py --real OUT   |   python -I -B s12_t6_contam_r2.py --fixture FIXDIR OUT      (OUT must exist and be empty)
import argparse, hashlib, itertools, json, os, sys, time, zipfile
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp

ROOT = _release_paths.expand('${REPO_ROOT}/')
POPS = ['Purkinje cell', 'Bergmann glia']
ALLPOPS = ['Astrocyte', 'Bergmann glia', 'Endothelial cell', 'GABAergic interneuron 1', 'GABAergic interneuron 2', 'GABAergic interneuron 3', 'Granule cell', 'Microglia', 'Oligodendrocyte', 'Oligodendrocyte progenitor cell', 'Pericyte', 'Purkinje cell', 'Unipolar brush cell']
AGES = [5, 18, 24, 30]
NPC, FLOOR = 30, 20
REAL_NFEAT, REAL_NCELLS, REAL_NUNITS = 1414, 280295, 312
THRESH = [2.0, 3.0, 4.0, 6.0]
PRIMARY = 4.0
TOL_PC = 1e-8
S6 = _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/')
S7 = _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/')
S6FILES = ['SHA256SUMS.txt', 'feature_registry.tsv', 'transform_keys.json', 'transform_parameters.npz', 'normalized_log10k_features.npz', 'cell_metadata.tsv.gz', 'unit_metadata.tsv']
PIN_REC = S7 + 'RECORD.json'
PIN_MAN = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv')
REAL_PINS = {
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): "7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/SHA256SUMS.txt'): "79cfae1a1bbf9f90c407425636600e31bbe1d9913f61698c3473b738212bcff0",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/cell_metadata.tsv.gz'): "ee6c4bca63f237825d7d80389e759263dbb45477e8535abb015b350f0595dbe6",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/feature_registry.tsv'): "d49b138d3a806db64fb37b6994a30034569cd530bfe3b6a82becf10f367f8c24",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/normalized_log10k_features.npz'): "082f1639be27db8791529bcaae7710b9cca4f8bf85b3f563ddd31de9455d9507",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/transform_keys.json'): "a9f42ea00582c9e9174c732cb3c9ee9090507bbcad53cd442d23c3d5f28889b8",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/transform_parameters.npz'): "0e4b000d17a09e33b9277f8e19e00eabf9e8dab319cbd71e179ad4b42cdb3b3b",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/unit_metadata.tsv'): "d2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/RECORD.json'): "594459e0e5a1b3a30f2f7b8e284968036aae3fafbaf02647c06871d9c5ec2a2a",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Astrocyte/nucleus_metadata.tsv.gz'): "f9cfcec06738e72ae2899c63b61930ec8878078efc44af68e660ded676745cae",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_PCs.npy'): "ed535ac9cdb36c51498155641172de92d9b6daa65e611ffa81a705b1a573ea13",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): "23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/pca_model.npz'): "34c1e8e2074b3c8247b09fac35d64acd6070e9f745b6604efdac3175f9642596",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Endothelial_cell/nucleus_metadata.tsv.gz'): "3792f519d08c5f2f4777e776eb6c50abf66d9ab19c465bb11542a04a5f37f3bb",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__GABAergic_interneuron_1/nucleus_metadata.tsv.gz'): "1fd038538f1e2b3baaf663199be054454fe211ba3a2c01b550c8fea1aaaad7b0",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__GABAergic_interneuron_2/nucleus_metadata.tsv.gz'): "e959562158c5038a2c2d52c4c57d122b7ed93bdbd4d78befe373bba9a170d9e8",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__GABAergic_interneuron_3/nucleus_metadata.tsv.gz'): "4016a094caf813c4423c4697990901dac1b3a1f05f17e0fa0ff4dd8385cff720",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Granule_cell/nucleus_metadata.tsv.gz'): "bedc0551f29c10a28bb629d1d31853c04e0304abbdb96571df4dda03d5b8e8f9",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Microglia/nucleus_metadata.tsv.gz'): "3c076729fac62813b7a447ba4168ee96efca66f28f186e32dc44614e71c074fa",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Oligodendrocyte/nucleus_metadata.tsv.gz'): "95bdeecd5dd9c9f0a00a0cfaf975c3b4dc5fc653123387bff60e4cfccc971a3b",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Oligodendrocyte_progenitor_cell/nucleus_metadata.tsv.gz'): "bc07159d6a00ad693801dda5f5131590069c5217e05f0429698301789bffb44d",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Pericyte/nucleus_metadata.tsv.gz'): "9177a525f64e3f01d5687d6d438dc3e7d5c2da397fb4c19be2e40909d70d3a9f",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'): "f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): "0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/pca_model.npz'): "eaf1468f99ed7935acf94a1e58b1a45f1fe1333ae6cc541c6afc550e4d128284",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Unipolar_brush_cell/nucleus_metadata.tsv.gz'): "60320fdf6fda6a286ab2df6b63a625a4f457ef0048d0c3a8c4abb6302ea29fd6",
}
RATIO_SCALE = 'linear CP10k (expm1 of the S06 log1p values)'
SCORE_SCALE = 'sum of S06 log1p(CP10k) values'
CHECKS = []


def stop(msg): raise SystemExit('STOP: ' + msg)


def ok(name, n=1): CHECKS.append((name, n))


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()


def exact_perm(wt, ki):
    if len(wt) < 2 or len(ki) < 2: return None
    v = np.array(list(wt) + list(ki), float); n = len(v); k = len(ki); obs = v[len(wt):].mean() - v[:len(wt)].mean(); st = []
    for idx in itertools.combinations(range(n), k):
        m = np.zeros(n, bool); m[list(idx)] = True; st.append(v[m].mean() - v[~m].mean())
    st = np.array(st); tol = 1e-12 * max(1.0, float(np.abs(st).max()))
    return dict(n_labelings=len(st), p_KI_higher=float((st >= obs - tol).mean()), p_KI_lower=float((st <= obs + tol).mean()), p_two_sided=float((np.abs(st) >= abs(obs) - tol).mean()))


def separation(wt, ki):
    if min(ki) > max(wt): return 'KI_higher'
    if max(ki) < min(wt): return 'KI_lower'
    return 'none'


def geno(s):
    k = str(s).strip().upper()
    if k in ('WT', 'KI'): return k
    stop('unrecognised genotype label: ' + repr(s))


def stem(name): return name.replace(' ', '_')


def p7(name, f): return S7 + 'discovery__%s/%s' % (stem(name), f)


def npz_shape(p):
    with zipfile.ZipFile(p) as z:
        with z.open('shape.npy') as f: return tuple(int(x) for x in np.lib.format.read_array(f))


def tkey(t): return ('%g' % t).replace('.', 'p')


def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True); g.add_argument('--real', action='store_true'); g.add_argument('--fixture'); ap.add_argument('out')
    a = ap.parse_args(); t0 = time.time(); real = a.real
    if real: base, pins = '', REAL_PINS
    else: base = a.fixture.rstrip('/') + '/'; pins = json.load(open(base + 'FIXTURE_PINS.json'))
    if not os.path.isdir(a.out) or os.listdir(a.out): stop('the output directory must exist and be empty')
    out = a.out.rstrip('/') + '/'
    for p, h in pins.items():
        if not os.path.isfile(base + p): stop('missing input ' + p)
        if sha256(base + p) != h: stop('sha256 differs from the pin: ' + p)
    keys = json.load(open(base + S6 + 'transform_keys.json')); pops = [kk['population'] for kk in keys if kk['mode'] == 'discovery']
    if len(set(pops)) != len(pops) or not all(n in pops for n in POPS): stop('transform_keys.json: a discovery population is listed twice, or Purkinje cell / Bergmann glia is missing')
    if real and sorted(pops) != sorted(ALLPOPS): stop('the discovery populations of transform_keys.json are not the 13 expected')
    want = [S6 + x for x in S6FILES] + [PIN_REC, PIN_MAN] + [p7(n, 'nucleus_metadata.tsv.gz') for n in pops] + [p7(n, f) for n in POPS for f in ('nucleus_PCs.npy', 'pca_model.npz')]
    if sorted(want) != sorted(pins): stop('the pin list is not exactly the %d inputs of this script' % len(want))
    ok('all %d inputs match their pinned sha256' % len(pins), len(pins))
    # ---------------- manifest
    man = pd.read_csv(base + PIN_MAN, sep='\t', dtype=str, keep_default_na=False)
    man = man[man['included'].str.strip().str.lower().isin(['true', '1', 'yes', 'y', 't'])].copy(); man['animal_id'] = man['animal_id'].str.strip(); man['genotype'] = man['genotype'].map(geno); man['age_weeks'] = man['age_weeks'].astype(float).astype(int)
    if len(man) != 24 or man['animal_id'].duplicated().any() or not (man.groupby(['age_weeks', 'genotype']).size() == 3).all(): stop('the pinned manifest is not 24 unique mice, 3 WT + 3 KI per age')
    man = man.set_index('animal_id'); ok('pinned manifest: 24 included mice, 3 WT + 3 KI per age', 24)
    # ---------------- proof (a): the S06 checksum chain
    sums = open(base + S6 + 'SHA256SUMS.txt').read().splitlines(); listed = {l.split(maxsplit=1)[1].strip(): l.split()[0] for l in sums if l.strip()}
    rec = json.load(open(base + PIN_REC))
    if sha256(base + S6 + 'SHA256SUMS.txt') != rec.get('source_transform_registry_sha256'): stop('S06 SHA256SUMS.txt does not equal source_transform_registry_sha256 of the S07 core RECORD.json')
    for f in S6FILES[1:]:
        if listed.get(f) != pins[S6 + f]: stop('S06 SHA256SUMS.txt does not list the pinned hash of ' + f)
    ok('proof (a): the S07 core RECORD pins the S06 stage A checksum list, which lists the pinned hashes of the registry, transform parameters, transform keys, normalized features, cell metadata and unit metadata', 7)
    # ---------------- proof (b), (c), (d): registry, units, cells, matrix shape
    reg = pd.read_csv(base + S6 + 'feature_registry.tsv', sep='\t', dtype={'var_position': int}); nf = len(reg)
    if (real and nf != REAL_NFEAT) or not (np.diff(reg['var_position'].to_numpy()) > 0).all() or reg['stable_id'].duplicated().any(): stop('proof (b): feature_registry.tsv does not have the expected row count with strictly increasing var_position and unique stable_id')
    ok('proof (b): feature_registry.tsv has %d rows, var_position strictly increasing, stable_id unique' % nf, nf)
    units = pd.read_csv(base + S6 + 'unit_metadata.tsv', sep='\t', dtype={'animal_id': str}); cells = pd.read_csv(base + S6 + 'cell_metadata.tsv.gz', sep='\t')
    for c in ('animal_id', 'cell_type_raw', 'age_weeks', 'genotype', 'n_nuclei', 'unit_index'):
        if c not in units.columns: stop('unit_metadata.tsv lacks column ' + c)
    for c in ('obs_name', 'unit_index', 'full_count_total'):
        if c not in cells.columns: stop('cell_metadata.tsv.gz lacks column ' + c)
    nU, N = len(units), len(cells)
    if not np.array_equal(units['unit_index'].to_numpy(), np.arange(nU)) or (real and nU != REAL_NUNITS): stop('proof (c): unit_metadata unit_index is not arange(len(units)) with the expected number of units')
    unit = cells['unit_index'].to_numpy()
    if unit.min() < 0 or unit.max() >= nU or cells['obs_name'].duplicated().any(): stop('proof (c): cell_metadata has a unit_index outside the units or a duplicated obs_name')
    if not np.array_equal(np.bincount(unit, minlength=nU), units['n_nuclei'].to_numpy()): stop('proof (c): bincount(cell_metadata.unit_index) differs from unit_metadata n_nuclei')
    shp = npz_shape(base + S6 + 'normalized_log10k_features.npz')
    if shp != (N, nf) or (real and shp != (REAL_NCELLS, REAL_NFEAT)): stop('proof (c): normalized_log10k_features.npz shape %s differs from (len(cell_metadata), len(feature_registry)) = (%d, %d)' % (shp, N, nf))
    ok('proof (c): unit_index is arange(%d); bincount(unit_index) equals n_nuclei; the normalized matrix shape %s equals (len(cell_metadata), len(feature_registry))' % (nU, shp), N)
    units['animal_id'] = units['animal_id'].str.strip(); units['g'] = units['genotype'].map(geno); units['age'] = units['age_weeks'].astype(float).astype(int)
    if not set(units['animal_id']) <= set(man.index) or len(set(units['animal_id'])) != 24: stop('proof (d): the units are not exactly the 24 included mice of the pinned manifest')
    for r_ in units.itertuples():
        if r_.g != man.loc[r_.animal_id, 'genotype'] or r_.age != man.loc[r_.animal_id, 'age_weeks']: stop('proof (d): unit_metadata disagrees with the pinned manifest for mouse ' + r_.animal_id)
    ok('proof (d): the units hold exactly the 24 included mice; genotype and age of every unit equal the pinned manifest', nU)
    X = sp.load_npz(base + S6 + 'normalized_log10k_features.npz').tocsr()
    if X.shape != (N, nf) or X.dtype != np.float64 or not np.isfinite(X.data).all(): stop('proof (c): the loaded normalized matrix differs from its stored shape, is not float64, or has a non-finite value')
    ok('normalized matrix loaded sparse (%d x %d, %d stored values, float64, all finite)' % (N, nf, X.nnz), X.nnz)
    # ---------------- proof (e): the S07 row mapping for every population
    pmask, rows = {}, {}
    for n in pops:
        pmask[n] = (units['cell_type_raw'].eq(n) & units['n_nuclei'].ge(FLOOR)).to_numpy(); rows[n] = np.flatnonzero(pmask[n][unit])
        sub = cells.iloc[rows[n]].reset_index(drop=True); nm = pd.read_csv(base + p7(n, 'nucleus_metadata.tsv.gz'), sep='\t')
        if list(nm.columns) != list(cells.columns) or len(nm) != len(sub) or nm['obs_name'].tolist() != sub['obs_name'].tolist() or not np.array_equal(nm['unit_index'].to_numpy(), sub['unit_index'].to_numpy()) or not np.array_equal(nm['full_count_total'].to_numpy(), sub['full_count_total'].to_numpy()): stop('proof (e): cells.iloc[all_cells] differs from the pinned nucleus_metadata of ' + n)
    allrows = np.concatenate([rows[n] for n in pops])
    if len(np.unique(allrows)) != len(allrows): stop('proof (e): a nucleus belongs to two populations')
    ok('proof (e): cells.iloc[flatnonzero(p_mask[unit])] equals the pinned nucleus_metadata of each of the %d populations (obs_name order, unit_index, full_count_total)' % len(pops), len(allrows))
    # ---------------- proof (f): the matrix rows of the scored populations reproduce the stored PC scores
    t_p = np.load(base + S6 + 'transform_parameters.npz'); elig = t_p['eligible']
    if elig.shape != (len(keys), nf): stop('proof (f): eligible has shape %s for %d keys and %d features' % (elig.shape, len(keys), nf))
    for name in POPS:
        ks = [i for i, kk in enumerate(keys) if kk['mode'] == 'discovery' and kk['population'] == name]
        pm = np.load(base + p7(name, 'pca_model.npz')); fi, rot, mean, scale = pm['feature_indices'], pm['rotation'], pm['mean'], pm['scale']; stored = np.load(base + p7(name, 'nucleus_PCs.npy'))
        if not np.array_equal(fi, np.flatnonzero(elig[ks[0]])) or rot.shape != (len(fi), NPC) or fi.min() < 0 or fi.max() >= nf: stop('proof (f): %s feature_indices differ from flatnonzero(eligible[key]) or rotation has the wrong shape' % name)
        if not np.array_equal(mean, t_p['mean'][ks[0]][fi]) or not np.array_equal(scale, t_p['scale'][ks[0]][fi]): stop('proof (f): %s mean or scale of pca_model differs from transform_parameters[key][eligible]' % name)
        if stored.shape != (len(rows[name]), NPC): stop('proof (f): %s nucleus_PCs has shape %s for %d nuclei' % (name, stored.shape, len(rows[name])))
        full = np.asarray(X[rows[name]][:, fi].multiply(1 / scale) @ rot) - (mean / scale) @ rot; err = np.abs(full - stored) / np.maximum(1.0, np.abs(stored))
        if not (err <= TOL_PC).all(): stop('proof (f): %s recomputed PC scores differ from nucleus_PCs by more than 1e-8 (relative to max(1, |stored|))' % name)
        ok('proof (f): %s the PC scores recomputed with the S07 operation order from the matrix rows equal nucleus_PCs for all %d PCs (max error %.3g)' % (name, NPC, float(err.max())), len(rows[name]))
    # ---------------- WT means of every population over all features (unit sums first, then divide)
    Gu = sp.csr_matrix((np.ones(N), (unit, np.arange(N))), shape=(nU, N)); US = np.asarray((Gu @ X).todense()); Xl = sp.csr_matrix((np.expm1(X.data), X.indices, X.indptr), shape=X.shape); USL = np.asarray((Gu @ Xl).todense()); del Xl; nu = units['n_nuclei'].to_numpy(float); unit_wt = (units['g'] == 'WT').to_numpy()
    Mp, Lp, Lb, den = np.zeros((len(pops), nf)), np.zeros((len(pops), nf)), np.zeros((len(pops), nf)), []
    for q, n in enumerate(pops):
        wu = np.flatnonzero(pmask[n] & unit_wt)
        if len(wu) == 0: stop('population %s has no WT unit with at least %d nuclei' % (n, FLOOR))
        Mp[q] = US[wu].sum(axis=0) / nu[wu].sum()
        Lp[q] = USL[wu].sum(axis=0) / nu[wu].sum(); Lb[q] = (USL[wu] / nu[wu][:, None]).mean(axis=0)
        den.append(dict(ratio_scale=RATIO_SCALE, population=n, n_units=int(pmask[n].sum()), n_nuclei=int(nu[pmask[n]].sum()), n_WT_units=len(wu), n_WT_nuclei=int(nu[wu].sum())))
    ok('WT means (linear CP10k pooled and mouse-balanced; log1p pooled) over %d features for each of %d populations' % (nf, len(pops)), len(pops))
    # ---------------- markers of every population
    def ratio_of(o, m): return np.where(m > 0, np.divide(o, m, out=np.zeros_like(o), where=m > 0), np.where(o > 0, np.inf, np.nan))
    mark, gene_rows = {}, []
    for q, n in enumerate(pops):
        own, own_b = Lp[q], Lb[q]; oth = np.delete(Lp, q, axis=0); oth_b = np.delete(Lb, q, axis=0); onames = [p for p in pops if p != n]; mo, mo_b = oth.max(axis=0), oth_b.max(axis=0); own_g, mo_g = Mp[q], np.delete(Mp, q, axis=0).max(axis=0)
        ratio, ratio_b = ratio_of(own, mo), ratio_of(own_b, mo_b); zl = np.where(mo > 0, 'defined', np.where(own > 0, 'max_other_zero_own_positive', 'own_zero'))
        mark[n] = {t: (own > 0) & (ratio >= t) for t in THRESH}
        d = dict(ratio_scale=RATIO_SCALE, population=n, feature_index=np.arange(nf), stable_id=reg['stable_id'].to_numpy(), dataset_symbol=reg['dataset_symbol'].to_numpy(), core_symbols=reg['core_symbols'].to_numpy(), own_WT_mean=own, max_other_WT_mean=mo, max_other_population=[onames[j] for j in oth.argmax(axis=0)], ratio=ratio, ratio_label=zl,
                 own_WT_mean_mouse_balanced=own_b, max_other_WT_mean_mouse_balanced=mo_b, ratio_mouse_balanced=ratio_b, own_WT_mean_log1p=own_g, max_other_WT_mean_log1p=mo_g, ratio_log1p=ratio_of(own_g, mo_g))
        for t in THRESH: d['marker_' + tkey(t)] = mark[n][t]
        d['marker_mouse_balanced_' + tkey(PRIMARY)] = (own_b > 0) & (ratio_b >= PRIMARY); gene_rows.append(pd.DataFrame(d))
    ok('markers defined for each of %d populations at %d thresholds' % (len(pops), len(THRESH)), len(pops) * len(THRESH))
    # ---------------- scores
    count_rows, list_rows, mean_rows, test_rows = [], [], [], []
    for P in POPS:
        rw = rows[P]; XP = X[rw].tocsc(); mice = units[pmask[P]].sort_values(['age', 'g', 'animal_id']); code = np.full(nU, -1); code[mice['unit_index'].to_numpy()] = np.arange(len(mice)); rc = code[unit[rw]]; cnt = np.bincount(rc, minlength=len(mice)).astype(float)
        for Q in [p for p in pops if p != P]:
            for t in THRESH:
                before = mark[Q][t]; used = before & ~mark[P][t]; idx = np.flatnonzero(used)
                count_rows.append(dict(ratio_scale=RATIO_SCALE, score_scale=SCORE_SCALE, scored_population=P, marker_population=Q, threshold=t, threshold_role='PRIMARY' if t == PRIMARY else 'SECONDARY', n_markers=int(before.sum()), n_removed_markers_of_scored_population=int((before & mark[P][t]).sum()), n_markers_used=len(idx), status='scored' if len(idx) else 'no_markers'))
                for j in idx: list_rows.append(dict(ratio_scale=RATIO_SCALE, scored_population=P, marker_population=Q, threshold=t, feature_index=int(j), stable_id=reg['stable_id'].iloc[j], dataset_symbol=reg['dataset_symbol'].iloc[j], core_symbols=reg['core_symbols'].iloc[j], ratio=float(gene_rows[pops.index(Q)]['ratio'].iloc[j]), own_WT_mean=float(Lp[pops.index(Q), j])))
                if not len(idx): continue
                sc = np.asarray(XP[:, idx].sum(axis=1)).ravel(); mm = np.bincount(rc, weights=sc, minlength=len(mice)) / cnt
                for i_, rr in enumerate(mice.itertuples()):
                    mean_rows.append(dict(score_scale=SCORE_SCALE, scored_population=P, marker_population=Q, threshold=t, threshold_role='PRIMARY' if t == PRIMARY else 'SECONDARY', animal_id=rr.animal_id, genotype=rr.g, age_weeks=rr.age, n_nuclei=int(cnt[i_]), n_markers_used=len(idx), mean_score=float(mm[i_])))
    MM = pd.DataFrame(mean_rows, columns=['score_scale', 'scored_population', 'marker_population', 'threshold', 'threshold_role', 'animal_id', 'genotype', 'age_weeks', 'n_nuclei', 'n_markers_used', 'mean_score'])
    for (P, Q, t, age), sub in MM.groupby(['scored_population', 'marker_population', 'threshold', 'age_weeks'], sort=False):
        w = sub[sub['genotype'] == 'WT']; ki = sub[sub['genotype'] == 'KI']; wv = w['mean_score'].to_numpy(float); kv = ki['mean_score'].to_numpy(float)
        rr = dict(score_scale=SCORE_SCALE, scored_population=P, marker_population=Q, threshold=t, threshold_role=sub['threshold_role'].iloc[0], age_weeks=age, n_markers_used=int(sub['n_markers_used'].iloc[0]), n_WT=len(w), n_KI=len(ki))
        if len(w) >= 1 and len(ki) >= 1:
            rr.update(WT_min=wv.min(), WT_max=wv.max(), WT_mean=wv.mean(), KI_min=kv.min(), KI_max=kv.max(), KI_mean=kv.mean(), diff_KI_minus_WT=kv.mean() - wv.mean(), separation=separation(wv, kv),
                      WT_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(w['animal_id'], wv)), KI_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(ki['animal_id'], kv)))
            pr = exact_perm(wv, kv)
            if pr: rr.update(pr)
        test_rows.append(rr)
    TT = pd.DataFrame(test_rows); srows = []
    if len(TT):
        for (P, Q, t), sub in TT.groupby(['scored_population', 'marker_population', 'threshold'], sort=False):
            row = dict(ratio_scale=RATIO_SCALE, score_scale=SCORE_SCALE, scored_population=P, marker_population=Q, threshold=t, threshold_role=sub['threshold_role'].iloc[0], n_markers_used=int(sub['n_markers_used'].iloc[0]))
            for r_ in sub.itertuples(): row.update({'diff_%dW' % r_.age_weeks: r_.diff_KI_minus_WT, 'separation_%dW' % r_.age_weeks: r_.separation, 'p_two_sided_%dW' % r_.age_weeks: r_.p_two_sided})
            srows.append(row)
    ok('contamination scores: %d (scored population, marker population, threshold) combinations scored, %d KI-vs-WT tests' % (len(MM.groupby(['scored_population', 'marker_population', 'threshold'])) if len(MM) else 0, len(TT)), len(TT))
    def wr(df, nm): df.to_csv(out + nm, sep='\t', index=False, float_format='%.10g', na_rep='NA', mode='x')
    wr(pd.concat(gene_rows, ignore_index=True), 'marker_genes.tsv'); wr(pd.DataFrame(count_rows), 'marker_counts.tsv'); wr(pd.DataFrame(list_rows, columns=['ratio_scale', 'scored_population', 'marker_population', 'threshold', 'feature_index', 'stable_id', 'dataset_symbol', 'core_symbols', 'ratio', 'own_WT_mean']), 'marker_lists.tsv'); wr(pd.DataFrame(den), 'wt_denominators.tsv')
    wr(MM, 'contamination_mouse_means.tsv'); wr(TT, 'contamination_kivswt.tsv'); wr(pd.DataFrame(srows), 'contamination_summary.tsv'); wr(pd.DataFrame(CHECKS, columns=['check_passed', 'items']), 'checks.tsv')
    CT = pd.DataFrame(count_rows)
    res = dict(mode='real' if real else 'fixture', checks=len(CHECKS), populations=pops, thresholds=THRESH, primary_threshold=PRIMARY, tests=len(TT), mouse_mean_rows=len(MM), marker_counts={('%s|%s|%g' % (r_.scored_population, r_.marker_population, r_.threshold)): [int(r_.n_markers), int(r_.n_removed_markers_of_scored_population), int(r_.n_markers_used)] for r_ in CT.itertuples()},
               n_markers_per_population={n: {('%g' % t): int(mark[n][t].sum()) for t in THRESH} for n in pops}, n_combinations_without_markers=int((CT['status'] == 'no_markers').sum()),
               ratio_scale=RATIO_SCALE, score_scale=SCORE_SCALE, definitions=dict(marker='own linear-scale WT mean > 0 and own / max over the other populations >= threshold (max_other 0 with own > 0: ratio inf); the pooled WT mean over all WT nuclei of the population, on the linear CP10k scale (expm1 of the S06 values), defines the markers; the log1p-scale ratio is written beside and unused', score='per nucleus the SUM of normalized expression of the used markers; per mouse the mean; KI vs WT within age, exact permutation test',
                                removal='genes that are markers of the scored population (same threshold) are removed from every Q list; the count is written', thresholds='2, 3, 4 (PRIMARY, a proposal), 6'),
               note='lists, counts and descriptive differences only; no interpretation; descriptive, not counted (S11 = C2)', python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, script_sha256=sha256(os.path.abspath(__file__)), seconds=round(time.time() - t0, 2))
    json.dump(res, open(out + 'RESULT.json', 'x'), indent=1, sort_keys=True); print(out)


if __name__ == '__main__': main()
#%% debug
