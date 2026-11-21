
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# S12 T5 , r1: EXACT GENE-SET DECOMPOSITION of the 12 discovery PC scores (Purkinje cell PC 5, 6, 11, 30; Bergmann glia PC 2, 4, 5, 6, 8, 9, 20, 29). DESCRIPTIVE, READ-ONLY, no interpretation, no refit.
#
# IDENTITY (s07_prepare_core_pca.py line 55): score[n, k] = sum over the population's features g of z[n, g] * rotation[g, k], z = X[n, g] / scale[g] - mean[g] / scale[g]. The sum is split by gene set; the parts add up to the full score.
# ROW MAPPING (proved here, STOP on failure): S07 takes the nuclei of a population as all_cells = flatnonzero(p_mask[unit]) with unit = cell_metadata.unit_index (s07_prepare_core_pca.py lines 39, 41, 47-48; the transform is line 55, pca_model.npz is written at line 57) and writes cells.iloc[all_cells] as nucleus_metadata (line 59);
#   S06 writes the matrix rows and cell_metadata rows in one loop order (s06_build.py lines 79-90), unit_index = arange(len(units)) (line 45), bincount(unit_index) == n_nuclei (line 87). Here: the matrix has len(cell_metadata) rows and len(feature_registry) columns, unit_index is arange,
#   bincount equals n_nuclei, cells.iloc[all_cells] equals the pinned nucleus_metadata of EVERY S07 population (obs_name order, unit_index, full_count_total), and the score recomputed with the S07 operation order equals the pinned nucleus_PCs for all 30 PCs: |recomputed - stored| <= 1e-8 * max(1, |stored|).
# POPULATION = the discovery keys of transform_keys.json (13); its nuclei = the S07 rule above (units with n_nuclei >= 20). WT nuclei = those whose unit is WT (unit_metadata genotype, equal to the pinned manifest for the 24 included mice).
# GENE SETS, per focus population P over its feature_indices (flatnonzero(eligible[key]), proved): own = WT mean of normalized expression in P; max_other = max over the other 12 populations of their WT means (all 13 populations' WT means are taken over the same 1,414 columns).
#   Two WT means are written: POOLED (all WT nuclei of the population pooled; DEFINES THE SETS, pending review's word) and MOUSE_BALANCED (mean over WT mice of the per-mouse mean; the weighting S07 uses for its fit). ratio = own / max_other; max_other == 0 and own > 0: ratio = +inf, label max_other_zero_own_positive; both 0: ratio = NaN, label both_zero (NOT intrinsic).
#   INTRINSIC(c) = ratio >= c. AMBIENT(P) = gene has potential_ambient_suspect TRUE or ambient_flagged_hit TRUE in ANY S11 release row of population P (labels must be exactly TRUE/FALSE, else STOP). Sets at cut-off c: INTRINSIC_CLEAN = INTRINSIC(c) and not AMBIENT; AMBIENT_FLAGGED = INTRINSIC(c) and AMBIENT (the gene is removed from INTRINSIC);
#   OTHER = everything else (a flagged gene that is not INTRINSIC stays in OTHER). Cut-offs 1 (PRIMARY), 1.5, 2, 3 are all decomposed in this run; cut-off 1 is the primary.
# PER MOUSE: the mean of every part (and of the full score) per PC; KI vs WT within age with the T3 exact test (exact_perm / separation copied from s12_t4_pc_r1.py); per-mouse full means must equal the T3 table; part means must add up to the full mean.
# Writes: gene_sets.tsv, wt_denominators.tsv, decomposition_mouse_means.tsv, decomposition_kivswt.tsv, decomposition_summary.tsv, checks.tsv, RESULT.json (nothing is overwritten; the output directory must be empty).
# usage: python -I -B s12_t5_decomp_r1.py --real OUT   |   python -I -B s12_t5_decomp_r1.py --fixture FIXDIR OUT      (OUT must exist and be empty)
import argparse, hashlib, itertools, json, os, sys, time, zipfile
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp

ROOT = _release_paths.expand('${REPO_ROOT}/')
POPS = [('Purkinje cell', [5, 6, 11, 30]), ('Bergmann glia', [2, 4, 5, 6, 8, 9, 20, 29])]       # (S07 / release registry name, PCs)
ALLPOPS = ['Astrocyte', 'Bergmann glia', 'Endothelial cell', 'GABAergic interneuron 1', 'GABAergic interneuron 2', 'GABAergic interneuron 3', 'Granule cell', 'Microglia', 'Oligodendrocyte', 'Oligodendrocyte progenitor cell', 'Pericyte', 'Purkinje cell', 'Unipolar brush cell']
AGES = [5, 18, 24, 30]
NPC, FLOOR = 30, 20
REAL_NFEAT, REAL_NCELLS, REAL_NUNITS = 1414, 280295, 312
CUTOFFS = [1.0, 1.5, 2.0, 3.0]
SETS = ['INTRINSIC_CLEAN', 'AMBIENT_FLAGGED', 'OTHER']
TOL_PC, TOL_SUM = 1e-8, 1e-9
S6 = _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/')
S7 = _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/')
S6FILES = ['SHA256SUMS.txt', 'feature_registry.tsv', 'transform_keys.json', 'transform_parameters.npz', 'normalized_log10k_features.npz', 'cell_metadata.tsv.gz', 'unit_metadata.tsv']
PIN_REC = S7 + 'RECORD.json'
PIN_MAN = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv')
PIN_T3 = _release_paths.expand('${REPO_ROOT}/results/s12_lanes/mouse_table.tsv')
PIN_REL = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz')
REAL_PINS = {
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz'): "40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298",
    _release_paths.expand('${REPO_ROOT}/results/s12_lanes/mouse_table.tsv'): "94cdfb6f7213dfee8d9ad01f30c0b6f99404335007a234cb2c4dc7b6d3938b00",
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


def label(series, what):
    lab = series.str.strip().str.upper(); bad = sorted(set(lab) - {'TRUE', 'FALSE'})
    if bad: stop('unexpected %s labels (not exactly TRUE/FALSE): %s' % (what, ', '.join(repr(b) for b in bad)))
    return (lab == 'TRUE').to_numpy()


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
    if len(set(pops)) != len(pops) or not all(n in pops for n, _ in POPS): stop('transform_keys.json: a discovery population is listed twice, or Purkinje cell / Bergmann glia is missing')
    if real and sorted(pops) != sorted(ALLPOPS): stop('the discovery populations of transform_keys.json are not the 13 expected')
    want = [S6 + x for x in S6FILES] + [PIN_REC, PIN_MAN, PIN_T3, PIN_REL] + [p7(n, 'nucleus_metadata.tsv.gz') for n in pops] + [p7(n, f) for n, _ in POPS for f in ('nucleus_PCs.npy', 'pca_model.npz')]
    if sorted(want) != sorted(pins): stop('the pin list is not exactly the %d inputs of this script' % len(want))
    ok('all %d inputs match their pinned sha256' % len(pins), len(pins))
    # ---------------- manifest, T3 table
    man = pd.read_csv(base + PIN_MAN, sep='\t', dtype=str, keep_default_na=False)
    man = man[man['included'].str.strip().str.lower().isin(['true', '1', 'yes', 'y', 't'])].copy(); man['animal_id'] = man['animal_id'].str.strip(); man['genotype'] = man['genotype'].map(geno); man['age_weeks'] = man['age_weeks'].astype(float).astype(int)
    if len(man) != 24 or man['animal_id'].duplicated().any() or not (man.groupby(['age_weeks', 'genotype']).size() == 3).all(): stop('the pinned manifest is not 24 unique mice, 3 WT + 3 KI per age')
    man = man.set_index('animal_id'); ok('pinned manifest: 24 included mice, 3 WT + 3 KI per age', 24)
    t3 = pd.read_csv(base + PIN_T3, sep='\t', dtype={'animal_id': str})
    # ---------------- proof (a): the S06 checksum chain
    sums = open(base + S6 + 'SHA256SUMS.txt').read().splitlines(); listed = {l.split(maxsplit=1)[1].strip(): l.split()[0] for l in sums if l.strip()}
    rec = json.load(open(base + PIN_REC))
    if sha256(base + S6 + 'SHA256SUMS.txt') != rec.get('source_transform_registry_sha256'): stop('S06 SHA256SUMS.txt does not equal source_transform_registry_sha256 of the S07 core RECORD.json')
    for f in S6FILES[1:]:
        if listed.get(f) != pins[S6 + f]: stop('S06 SHA256SUMS.txt does not list the pinned hash of ' + f)
    ok('proof (a): the S07 core RECORD pins the S06 stage A checksum list, which lists the pinned hashes of the registry, transform parameters, transform keys, normalized features, cell metadata and unit metadata', 7)
    # ---------------- proof (b): registry, units, cells, matrix shape
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
    # ---------------- WT means of every population over all features (unit sums first, then divide: exact for equal values)
    Gu = sp.csr_matrix((np.ones(N), (unit, np.arange(N))), shape=(nU, N)); US = np.asarray((Gu @ X).todense()); nu = units['n_nuclei'].to_numpy(float); unit_wt = (units['g'] == 'WT').to_numpy()
    Mp, Mb, den = np.zeros((len(pops), nf)), np.zeros((len(pops), nf)), []
    for q, n in enumerate(pops):
        wu = np.flatnonzero(pmask[n] & unit_wt)
        if len(wu) == 0: stop('population %s has no WT unit with at least %d nuclei' % (n, FLOOR))
        Mp[q] = US[wu].sum(axis=0) / nu[wu].sum(); Mb[q] = (US[wu] / nu[wu][:, None]).mean(axis=0)
        den.append(dict(population=n, n_units=int(pmask[n].sum()), n_nuclei=int(nu[pmask[n]].sum()), n_WT_units=len(wu), n_WT_nuclei=int(nu[wu].sum())))
    ok('WT means (pooled and mouse-balanced) over %d features for each of %d populations' % (nf, len(pops)), len(pops))
    # ---------------- release flags
    rel = pd.read_csv(base + PIN_REL, sep='\t', dtype=str, keep_default_na=False, usecols=['stable_id', 'population_registry_name', 'potential_ambient_suspect', 'ambient_flagged_hit'])
    rel['sus'] = label(rel['potential_ambient_suspect'], 'potential_ambient_suspect'); rel['hit'] = label(rel['ambient_flagged_hit'], 'ambient_flagged_hit')
    t_p = np.load(base + S6 + 'transform_parameters.npz'); elig = t_p['eligible']
    if elig.shape != (len(keys), nf): stop('proof (f): eligible has shape %s for %d keys and %d features' % (elig.shape, len(keys), nf))
    gene_rows, mean_rows, test_rows, res_pop = [], [], [], {}
    for name, pcs_sel in POPS:
        iP = pops.index(name); ks = [i for i, kk in enumerate(keys) if kk['mode'] == 'discovery' and kk['population'] == name]
        pm = np.load(base + p7(name, 'pca_model.npz')); fi, rot, mean, scale = pm['feature_indices'], pm['rotation'], pm['mean'], pm['scale']; stored = np.load(base + p7(name, 'nucleus_PCs.npy'))
        if not np.array_equal(fi, np.flatnonzero(elig[ks[0]])) or rot.shape != (len(fi), NPC) or fi.min() < 0 or fi.max() >= nf: stop('proof (f): %s feature_indices differ from flatnonzero(eligible[key]) or rotation has the wrong shape' % name)
        if not np.array_equal(mean, t_p['mean'][ks[0]][fi]) or not np.array_equal(scale, t_p['scale'][ks[0]][fi]): stop('proof (f): %s mean or scale of pca_model differs from transform_parameters[key][eligible]' % name)
        rw = rows[name]
        if stored.shape != (len(rw), NPC): stop('proof (f): %s nucleus_PCs has shape %s for %d nuclei' % (name, stored.shape, len(rw)))
        Xsc = sp.csc_matrix(X[rw][:, fi].multiply(1 / scale)); base_ = (mean / scale) @ rot
        full = np.asarray(Xsc @ rot) - base_
        err = np.abs(full - stored) / np.maximum(1.0, np.abs(stored))
        if not (err <= TOL_PC).all(): stop('proof (f): %s recomputed PC scores differ from nucleus_PCs by more than 1e-8 (relative to max(1, |stored|))' % name)
        ok('proof (f): %s feature_indices equal flatnonzero(eligible[%d]); mean and scale equal transform_parameters; the PC scores recomputed with the S07 operation order equal nucleus_PCs for all %d PCs (max error %.3g)' % (name, ks[0], NPC, float(err.max())), len(rw))
        # gene sets
        own, own_b = Mp[iP, fi], Mb[iP, fi]; oth = np.delete(Mp, iP, axis=0)[:, fi]; oth_b = np.delete(Mb, iP, axis=0)[:, fi]; onames = [p for p in pops if p != name]
        mo, mo_b = oth.max(axis=0), oth_b.max(axis=0)
        def ratio_of(o, m): return np.where(m > 0, np.divide(o, m, out=np.zeros_like(o), where=m > 0), np.where(o > 0, np.inf, np.nan))
        ratio, ratio_b = ratio_of(own, mo), ratio_of(own_b, mo_b)
        zl = np.where(mo > 0, 'defined', np.where(own > 0, 'max_other_zero_own_positive', 'both_zero'))
        R_ = rel[rel['population_registry_name'] == name]
        if R_.empty: stop(name + ': no rows of this population in the S11 release')
        sid = reg['stable_id'].to_numpy()[fi]; sus_ids, hit_ids, in_ids = set(R_.loc[R_['sus'], 'stable_id']), set(R_.loc[R_['hit'], 'stable_id']), set(R_['stable_id'])
        nrel = R_.groupby('stable_id').size(); sus = np.array([s in sus_ids for s in sid]); hit = np.array([s in hit_ids for s in sid]); amb = sus | hit
        def set_of(rt, c): return np.where(rt >= c, np.where(amb, 'AMBIENT_FLAGGED', 'INTRINSIC_CLEAN'), 'OTHER')
        sets = {c: set_of(ratio, c) for c in CUTOFFS}
        G = pd.DataFrame(dict(population=name, feature_index=fi, stable_id=sid, dataset_symbol=reg['dataset_symbol'].to_numpy()[fi], core_symbols=reg['core_symbols'].to_numpy()[fi], own_WT_mean=own, max_other_WT_mean=mo, max_other_population=[onames[j] for j in oth.argmax(axis=0)], ratio=ratio, zero_denominator=zl,
                              own_WT_mean_mouse_balanced=own_b, max_other_WT_mean_mouse_balanced=mo_b, max_other_population_mouse_balanced=[onames[j] for j in oth_b.argmax(axis=0)], ratio_mouse_balanced=ratio_b, in_release=[s in in_ids for s in sid], release_rows=[int(nrel.get(s, 0)) for s in sid],
                              potential_ambient_suspect_any=sus, ambient_flagged_hit_any=hit, set=sets[1.0], set_cutoff_1p5=sets[1.5], set_cutoff_2=sets[2.0], set_cutoff_3=sets[3.0], set_mouse_balanced_cutoff_1=set_of(ratio_b, 1.0)))
        gene_rows.append(G)
        # decomposition
        parts = {}
        for c in CUTOFFS:
            tot = np.zeros_like(full)
            for s in SETS:
                idx = np.flatnonzero(sets[c] == s)
                pt = (np.asarray(Xsc[:, idx] @ rot[idx]) - (mean / scale)[idx] @ rot[idx]) if len(idx) else np.zeros_like(full)
                parts[(c, s)] = pt; tot = tot + pt
            if np.abs(tot - full).max() > TOL_SUM: stop('%s cut-off %g: the parts do not add up to the full score (1e-9)' % (name, c))
        ok('%s: for each of %d cut-offs the three parts add up to the recomputed full score for all %d nuclei and %d PCs (absolute 1e-9)' % (name, len(CUTOFFS), len(rw), NPC), len(rw) * NPC * len(CUTOFFS))
        mice = units[pmask[name]].sort_values(['age', 'g', 'animal_id']); code = np.full(nU, -1); code[mice['unit_index'].to_numpy()] = np.arange(len(mice)); rc = code[unit[rw]]; cnt = np.bincount(rc, minlength=len(mice)).astype(float)
        T = t3[t3['population'] == stem(name)].set_index('animal_id')
        if set(T.index) != set(mice['animal_id']): stop(name + ': the T3 table and the units hold different mice')
        def mmean(v): return np.bincount(rc, weights=v, minlength=len(mice)) / cnt
        for k in pcs_sel:
            fm = mmean(full[:, k - 1])
            for i_, aid in enumerate(mice['animal_id']):
                t3v = float(T.loc[aid, 'discovery_PC%d_mean' % k])
                if abs(t3v - fm[i_]) > 1e-8 * max(1.0, abs(t3v)): stop('%s %s PC%d: the full PC mean differs from the T3 table value' % (name, aid, k))
            for c in CUTOFFS:
                if np.abs(sum(mmean(parts[(c, s)][:, k - 1]) for s in SETS) - fm).max() > TOL_SUM: stop('%s PC%d cut-off %g: the part means do not add up to the full mean' % (name, k, c))
            for part, c, vec in [('FULL', None, fm)] + [(s, c, mmean(parts[(c, s)][:, k - 1])) for c in CUTOFFS for s in SETS]:
                for i_, (aid, rr) in enumerate(zip(mice['animal_id'], mice.itertuples())):
                    mean_rows.append(dict(population=name, PC=k, part=part, cutoff='NA' if c is None else c, cutoff_role='NA' if c is None else ('PRIMARY' if c == 1.0 else 'SECONDARY'), animal_id=aid, genotype=rr.g, age_weeks=rr.age, n_nuclei=int(cnt[i_]), mean=float(vec[i_])))
        ok('%s: per-mouse full means equal the T3 table for %d PCs; part means add up to the full mean' % (name, len(pcs_sel)), len(mice) * len(pcs_sel))
        res_pop[name] = dict(n_features=len(fi), n_nuclei=len(rw), n_mice=len(mice), max_pc_error=float(err.max()), zero_denominator={l: int((zl == l).sum()) for l in ('defined', 'max_other_zero_own_positive', 'both_zero')}, n_suspect=int(sus.sum()), n_hit_flag=int(hit.sum()), n_flagged=int(amb.sum()), n_not_in_release=int((~np.array([s in in_ids for s in sid])).sum()),
                             sets_by_cutoff={('%g' % c): {s: int((sets[c] == s).sum()) for s in SETS} for c in CUTOFFS}, n_sets_differ_pooled_vs_balanced_cutoff_1=int((set_of(ratio_b, 1.0) != sets[1.0]).sum()))
    MM = pd.DataFrame(mean_rows)
    for name, pcs_sel in POPS:
        for k in pcs_sel:
            for age in AGES:
                for (part, c), sub in MM[(MM['population'] == name) & (MM['PC'] == k) & (MM['age_weeks'] == age)].groupby(['part', 'cutoff'], sort=False):
                    w = sub[sub['genotype'] == 'WT']; ki = sub[sub['genotype'] == 'KI']; wv = w['mean'].to_numpy(float); kv = ki['mean'].to_numpy(float)
                    rr = dict(population=name, PC=k, age_weeks=age, part=part, cutoff=c, cutoff_role=sub['cutoff_role'].iloc[0], n_WT=len(w), n_KI=len(ki))
                    if len(w) >= 1 and len(ki) >= 1:
                        rr.update(WT_min=wv.min(), WT_max=wv.max(), WT_mean=wv.mean(), KI_min=kv.min(), KI_max=kv.max(), KI_mean=kv.mean(), diff_KI_minus_WT=kv.mean() - wv.mean(), separation=separation(wv, kv),
                                  WT_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(w['animal_id'], wv)), KI_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(ki['animal_id'], kv)))
                        pr = exact_perm(wv, kv)
                        if pr: rr.update(pr)
                    test_rows.append(rr)
    TT = pd.DataFrame(test_rows); look = {(r_.population, r_.PC, r_.age_weeks, r_.part, str(r_.cutoff)): r_ for r_ in TT.itertuples()}; srows = []
    for name, pcs_sel in POPS:
        for k in pcs_sel:
            for age in AGES:
                for c in CUTOFFS:
                    row = dict(population=name, PC=k, age_weeks=age, cutoff=c, cutoff_role='PRIMARY' if c == 1.0 else 'SECONDARY'); f_ = look[(name, k, age, 'FULL', 'NA')]; row.update(n_WT=f_.n_WT, n_KI=f_.n_KI)
                    for part in ['FULL'] + SETS:
                        t_ = look[(name, k, age, part, 'NA' if part == 'FULL' else str(c))]; pre = part
                        row.update({pre + '_diff': getattr(t_, 'diff_KI_minus_WT', np.nan), pre + '_separation': getattr(t_, 'separation', ''), pre + '_p_two_sided': getattr(t_, 'p_two_sided', np.nan)})
                    for s in SETS: row['n_genes_' + s] = res_pop[name]['sets_by_cutoff']['%g' % c][s]
                    srows.append(row)
    def wr(df, nm): df.to_csv(out + nm, sep='\t', index=False, float_format='%.10g', na_rep='NA', mode='x')
    wr(pd.concat(gene_rows, ignore_index=True), 'gene_sets.tsv'); wr(pd.DataFrame(den), 'wt_denominators.tsv'); wr(MM, 'decomposition_mouse_means.tsv'); wr(TT, 'decomposition_kivswt.tsv'); wr(pd.DataFrame(srows), 'decomposition_summary.tsv'); wr(pd.DataFrame(CHECKS, columns=['check_passed', 'items']), 'checks.tsv')
    res = dict(mode='real' if real else 'fixture', checks=len(CHECKS), populations=pops, focus=res_pop, cutoffs=CUTOFFS, primary_cutoff=1.0, tests=len(TT), mouse_mean_rows=len(MM),
               definitions=dict(wt_mean='pooled over all WT nuclei of the population (units with n_nuclei >= 20); the mouse-balanced mean is written beside it; the sets use the pooled mean', ratio='own / max over the other populations; max_other 0 and own > 0 = inf; both 0 = NaN (not intrinsic)', intrinsic='ratio >= cut-off',
                                ambient='potential_ambient_suspect TRUE or ambient_flagged_hit TRUE in any S11 release row of the population', sets='INTRINSIC_CLEAN = intrinsic and not ambient; AMBIENT_FLAGGED = intrinsic and ambient; OTHER = the rest'),
               release_labels_checked=['potential_ambient_suspect', 'ambient_flagged_hit'], note='lists, counts and descriptive differences only; no interpretation; descriptive, not counted (S11 = C2)', python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
               script_sha256=sha256(os.path.abspath(__file__)), seconds=round(time.time() - t0, 2))
    json.dump(res, open(out + 'RESULT.json', 'x'), indent=1, sort_keys=True); print(out)


if __name__ == '__main__': main()
#%% debug
