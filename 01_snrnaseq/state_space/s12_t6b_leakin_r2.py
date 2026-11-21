
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# S12 T6b leak-in sensitivity , r2: full-source-gene marker universe; code and synthetic fixtures only until review byte review and review.
# READ-ONLY/descriptive. Marker rule is declared before any real-data execution: pooled WT mean >=1 linear CP10k, detection in >=20% of that type's WT nuclei, ratio >=2/3/4/6 against every other population, top 50 per population/threshold by pooled ratio, feature-index tie-break.
# Granule and GABAergic interneuron 1 marker-set measurability is explicitly reported at every threshold. Scores are per mouse; KI-vs-WT comparisons remain within age.
# Synthetic harness MUST replace ROOT, REAL_PINS and REAL_CFG immediately after exec and call main only with --fixture paths. No real-data execution before review.
# usage: python -I -B s12_t6b_leakin_r2.py --real OUT | --fixture FIXDIR OUT (OUT must exist and be empty)
import argparse
import hashlib
import itertools
import json
import os
import sys
import time
import zipfile

import h5py
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
from anndata.io import read_elem

ROOT = _release_paths.expand('${DATA_ROOT}/')
H5AD_REL = 'data/adata_leiden_0.175_annotated_cell_types_v3.h5ad'
S6 = 'results/s06_stageA_v01/'
S7 = 'results/s07_core_representation_v01/'
PIN_REC = S7 + 'RECORD.json'
PIN_MAN = 'results/sca1_audit_run01/animal_manifest.tsv'
POPS = ['Purkinje cell', 'Bergmann glia']
ALLPOPS = ['Astrocyte', 'Bergmann glia', 'Endothelial cell', 'GABAergic interneuron 1', 'GABAergic interneuron 2', 'GABAergic interneuron 3', 'Granule cell', 'Microglia', 'Oligodendrocyte', 'Oligodendrocyte progenitor cell', 'Pericyte', 'Purkinje cell', 'Unipolar brush cell']
MEASURABILITY_TARGETS = ['Granule cell', 'GABAergic interneuron 1']
AGES = [5, 18, 24, 30]
FLOOR = 20
REAL_NOBS, REAL_NVAR, REAL_NCELLS, REAL_NUNITS = 318312, 26374, 280295, 312
THRESH = [2.0, 3.0, 4.0, 6.0]
PRIMARY = 4.0
MEAN_FLOOR_CP10K = 1.0
PREVALENCE_FLOOR = 0.20
CAP = 50
H5AD_SHA256 = '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'
S6_REQUIRED = ['feature_registry.tsv', 'cell_metadata.tsv.gz', 'unit_metadata.tsv']
P7_METADATA = [S7 + 'discovery__%s/nucleus_metadata.tsv.gz' % n.replace(' ', '_') for n in ALLPOPS]

# S06 hashes are retained to bind its panel metadata and row/unit maps to the pinned S07 record.
REAL_PINS = {
    H5AD_REL: H5AD_SHA256,
    'results/sca1_audit_run01/animal_manifest.tsv': '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07',
    S6 + 'SHA256SUMS.txt': '79cfae1a1bbf9f90c407425636600e31bbe1d9913f61698c3473b738212bcff0',
    S6 + 'feature_registry.tsv': 'd49b138d3a806db64fb37b6994a30034569cd530bfe3b6a82becf10f367f8c24',
    S6 + 'cell_metadata.tsv.gz': 'ee6c4bca63f237825d7d80389e759263dbb45477e8535abb015b350f0595dbe6',
    S6 + 'unit_metadata.tsv': 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657',
    PIN_REC: '594459e0e5a1b3a30f2f7b8e284968036aae3fafbaf02647c06871d9c5ec2a2a',
    P7_METADATA[0]: 'f9cfcec06738e72ae2899c63b61930ec8878078efc44af68e660ded676745cae',
    P7_METADATA[1]: '23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540',
    P7_METADATA[2]: '3792f519d08c5f2f4777e776eb6c50abf66d9ab19c465bb11542a04a5f37f3bb',
    P7_METADATA[3]: '1fd038538f1e2b3baaf663199be054454fe211ba3a2c01b550c8fea1aaaad7b0',
    P7_METADATA[4]: 'e959562158c5038a2c2d52c4c57d122b7ed93bdbd4d78befe373bba9a170d9e8',
    P7_METADATA[5]: '4016a094caf813c4423c4697990901dac1b3a1f05f17e0fa0ff4dd8385cff720',
    P7_METADATA[6]: 'bedc0551f29c10a28bb629d1d31853c04e0304abbdb96571df4dda03d5b8e8f9',
    P7_METADATA[7]: '3c076729fac62813b7a447ba4168ee96efca66f28f186e32dc44614e71c074fa',
    P7_METADATA[8]: '95bdeecd5dd9c9f0a00a0cfaf975c3b4dc5fc653123387bff60e4cfccc971a3b',
    P7_METADATA[9]: 'bc07159d6a00ad693801dda5f5131590069c5217e05f0429698301789bffb44d',
    P7_METADATA[10]: '9177a525f64e3f01d5687d6d438dc3e7d5c2da397fb4c19be2e40909d70d3a9f',
    P7_METADATA[11]: '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721',
    P7_METADATA[12]: '60320fdf6fda6a286ab2df6b63a625a4f457ef0048d0c3a8c4abb6302ea29fd6',
}
REAL_CFG = {'root': ROOT, 'pins': REAL_PINS}
RATIO_SCALE = 'raw counts normalized per nucleus to linear CP10k using the full-library count total'
SCORE_SCALE = 'per-nucleus sum of log1p(CP10k) over selected full-source-gene markers'
CHECKS = []


def stop(msg):
    raise SystemExit('STOP: ' + msg)


def ok(name, n=1):
    CHECKS.append((name, int(n)))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def exact_perm(wt, ki):
    if len(wt) != 3 or len(ki) != 3:
        return None
    values = np.asarray(list(wt) + list(ki), dtype=float)
    observed = values[3:].mean() - values[:3].mean()
    stats = []
    for selected in itertools.combinations(range(6), 3):
        mask = np.zeros(6, dtype=bool)
        mask[list(selected)] = True
        stats.append(values[mask].mean() - values[~mask].mean())
    stats = np.asarray(stats)
    tol = 1e-12 * max(1.0, float(np.abs(stats).max()))
    return {
        'n_labelings': int(len(stats)),
        'p_KI_higher': float((stats >= observed - tol).mean()),
        'p_KI_lower': float((stats <= observed + tol).mean()),
        'p_two_sided': float((np.abs(stats) >= abs(observed) - tol).mean()),
    }


def tkey(threshold):
    return ('%g' % threshold).replace('.', 'p')


def stem(name):
    return name.replace(' ', '_')


def p7(name):
    return S7 + 'discovery__%s/nucleus_metadata.tsv.gz' % stem(name)


def npz_shape(path):
    with zipfile.ZipFile(path) as archive:
        with archive.open('shape.npy') as stream:
            return tuple(int(x) for x in np.lib.format.read_array(stream))


def validate_manifest(base, units):
    manifest = pd.read_csv(base + PIN_MAN, sep='\t', dtype=str, keep_default_na=False)
    manifest = manifest[manifest['included'].str.strip().str.lower().isin(['true', '1', 'yes', 'y', 't'])].copy()
    manifest['animal_id'] = manifest['animal_id'].str.strip()
    labels = manifest['genotype'].astype(str).str.strip().str.upper()
    if not labels.isin(['WT', 'KI']).all():
        stop('unrecognised genotype label in the pinned manifest')
    manifest['genotype'] = labels
    manifest['age_weeks'] = manifest['age_weeks'].astype(float).astype(int)
    if len(manifest) != 24 or manifest['animal_id'].duplicated().any():
        stop('the pinned manifest is not 24 unique mice')
    if set(manifest['age_weeks']) != set(AGES) or not (manifest.groupby(['age_weeks', 'genotype']).size() == 3).all():
        stop('the pinned manifest does not contain 3 WT + 3 KI at each declared age')
    manifest = manifest.set_index('animal_id')
    units = units.copy()
    units['animal_id'] = units['animal_id'].astype(str).str.strip()
    units['genotype_norm'] = units['genotype'].astype(str).str.strip().str.upper()
    units['age_norm'] = units['age_weeks'].astype(float).astype(int)
    if not set(units['animal_id']) <= set(manifest.index) or len(set(units['animal_id'])) != 24:
        stop('unit_metadata does not refer exactly to the 24 included mice')
    for row in units.itertuples():
        if row.genotype_norm != manifest.loc[row.animal_id, 'genotype'] or row.age_norm != manifest.loc[row.animal_id, 'age_weeks']:
            stop('unit_metadata genotype/age differs from the pinned manifest for ' + row.animal_id)
    return manifest, units


def verify_inputs(base, pins):
    for rel, expected in pins.items():
        path = base + rel
        if not os.path.isfile(path):
            stop('missing pinned input ' + rel)
        if sha256(path) != expected:
            stop('sha256 differs from the pin: ' + rel)
    if len(pins) != len(set(pins)):
        stop('the pin list has duplicate paths')
    ok('all %d declared inputs match their pinned sha256' % len(pins), len(pins))


def load_inputs(base, real):
    cells = pd.read_csv(base + S6 + 'cell_metadata.tsv.gz', sep='\t')
    units = pd.read_csv(base + S6 + 'unit_metadata.tsv', sep='\t', dtype={'animal_id': str})
    for column in ('animal_id', 'cell_type_raw', 'age_weeks', 'genotype', 'n_nuclei', 'unit_index'):
        if column not in units.columns:
            stop('unit_metadata.tsv lacks column ' + column)
    for column in ('obs_name', 'unit_index', 'full_count_total'):
        if column not in cells.columns:
            stop('cell_metadata.tsv.gz lacks column ' + column)
    if not np.array_equal(units['unit_index'].to_numpy(), np.arange(len(units))):
        stop('unit_metadata unit_index is not arange(len(units))')
    if real and (len(units) != REAL_NUNITS or len(cells) != REAL_NCELLS):
        stop('pinned S06 cell/unit row counts differ from the declared real dimensions')
    if cells['obs_name'].duplicated().any():
        stop('cell_metadata.tsv.gz has duplicate obs_name')
    unit_index = cells['unit_index'].to_numpy(dtype=np.int64)
    if len(unit_index) and (unit_index.min() < 0 or unit_index.max() >= len(units)):
        stop('cell_metadata.unit_index is outside unit_metadata')
    if not np.array_equal(np.bincount(unit_index, minlength=len(units)), units['n_nuclei'].to_numpy(dtype=np.int64)):
        stop('bincount(cell_metadata.unit_index) differs from unit_metadata n_nuclei')
    sums_path = base + S6 + 'SHA256SUMS.txt'
    sums_lines = [line.split(maxsplit=1) for line in open(sums_path, encoding='utf-8') if line.strip()]
    listed = {parts[1].strip(): parts[0] for parts in sums_lines}
    record = json.load(open(base + PIN_REC, encoding='utf-8'))
    if sha256(sums_path) != record.get('source_transform_registry_sha256'):
        stop('S06 SHA256SUMS.txt does not equal source_transform_registry_sha256 of the S07 core RECORD.json')
    for rel in S6_REQUIRED:
        if listed.get(rel) != sha256(base + S6 + rel):
            stop('S06 SHA256SUMS.txt does not list the pinned hash of ' + rel)
    ok('S06-to-S07 checksum chain binds the panel registry and cell/unit metadata', len(S6_REQUIRED) + 1)
    manifest, units = validate_manifest(base, units)
    ok('pinned manifest: 24 included mice, 3 WT + 3 KI within each of four ages', 24)
    unit_index = cells['unit_index'].to_numpy(dtype=np.int64)
    unit_pop = units['cell_type_raw'].astype(str).to_numpy()
    unit_genotype = units['genotype_norm'].to_numpy()
    unit_age = units['age_norm'].to_numpy(dtype=np.int64)
    unit_n = units['n_nuclei'].to_numpy(dtype=np.int64)
    unit_active = unit_n >= FLOOR
    if not set(ALLPOPS) <= set(unit_pop):
        stop('unit_metadata does not contain all 13 declared cell populations')
    for name in ALLPOPS:
        wu = unit_active & (unit_pop == name) & (unit_genotype == 'WT')
        if int(wu.sum()) == 0:
            stop('population %s has no WT unit with at least %d nuclei' % (name, FLOOR))
    for name in POPS:
        for age in AGES:
            for genotype in ('WT', 'KI'):
                selected = unit_active & (unit_pop == name) & (unit_age == age) & (unit_genotype == genotype)
                if int(selected.sum()) != 3:
                    stop('%s age %d requires exactly 3 active %s mice for its per-mouse comparison' % (name, age, genotype))
    ok('all 13 populations have active WT units; both scored populations have 3 WT + 3 KI per age', len(ALLPOPS) + 2 * len(AGES) * 2)
    return cells, units, unit_index, unit_pop, unit_genotype, unit_age, unit_n, unit_active


def load_gene_registry(h5ad, panel_path):
    with h5py.File(h5ad, 'r') as handle:
        if 'obs' not in handle or 'var' not in handle or 'layers/counts' not in handle:
            stop('source H5AD must contain obs, var and layers/counts')
        obs = read_elem(handle['obs'])
        var = read_elem(handle['var'])
        matrix = handle['layers/counts']
        if tuple(matrix.attrs.get('shape', ())) != (len(obs), len(var)):
            stop('full-source counts shape differs from H5AD obs/var dimensions')
        if matrix.attrs.get('encoding-type') != 'csr_matrix':
            stop('full-source layers/counts must be a CSR matrix')
        if 'gene_ids' not in var.columns:
            stop('full-source var metadata lacks gene_ids')
        stable = var['gene_ids'].astype(str).to_numpy()
        symbols = var.index.astype(str).to_numpy()
    if len(stable) == 0 or pd.isna(stable).any() or (stable == '').any() or pd.Index(stable).duplicated().any():
        stop('full-source gene_ids are empty or duplicated')
    registry = pd.DataFrame({'feature_index': np.arange(len(stable), dtype=np.int64), 'stable_id': stable, 'dataset_symbol': symbols})
    panel = pd.read_csv(panel_path, sep='\t', dtype={'stable_id': str})
    if not {'stable_id', 'var_position'} <= set(panel.columns):
        stop('S06 panel feature_registry lacks stable_id or var_position')
    pos = panel['var_position'].to_numpy(dtype=np.int64)
    if len(panel) != 1414 or len(np.unique(pos)) != len(pos) or (pos < 0).any() or (pos >= len(registry)).any():
        stop('S06 panel registry is not the expected 1,414 unique positions in the all-gene source')
    if not np.array_equal(registry['stable_id'].to_numpy()[pos], panel['stable_id'].astype(str).to_numpy()):
        stop('S06 panel stable IDs do not map to their positions in the full-source gene registry')
    return registry


def source_row_map(h5ad, cells, units, real):
    with h5py.File(h5ad, 'r') as handle:
        obs = read_elem(handle['obs'])
        counts = handle['layers/counts']
        shape = tuple(int(x) for x in counts.attrs['shape'])
        if real and shape != (REAL_NOBS, REAL_NVAR):
            stop('real full-source count matrix dimensions differ from the pinned expected dimensions')
        row_index = obs.index.get_indexer(cells['obs_name'].astype(str))
        if (row_index < 0).any() or len(np.unique(row_index)) != len(cells) or (np.diff(row_index) <= 0).any():
            stop('S06 retained cell_metadata rows do not map uniquely and in source order to H5AD obs')
    source_to_cell = np.full(shape[0], -1, dtype=np.int64)
    source_to_cell[row_index] = np.arange(len(row_index), dtype=np.int64)
    return source_to_cell, shape


def iter_count_chunks(handle, shape, source_to_cell, cells, cell_unit, units, unit_pop, unit_genotype, unit_n, unit_active, linear_sums, detections, wt_nuclei, target_sums=None, target_nuclei=None, score_matrices=None, target_unit_maps=None):
    matrix = handle['layers/counts']
    data_ds, index_ds, ptr_ds = matrix['data'], matrix['indices'], matrix['indptr']
    n_source, n_features = shape
    chunk_rows = 2000
    full_totals = cells['full_count_total'].to_numpy(dtype=np.float64)
    for start in range(0, n_source, chunk_rows):
        end = min(n_source, start + chunk_rows)
        ptr = ptr_ds[start:end + 1].astype(np.int64)
        left, right = int(ptr[0]), int(ptr[-1])
        local_ptr = ptr - left
        values = data_ds[left:right].astype(np.float64)
        cols = index_ds[left:right].astype(np.int64)
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            stop('full-source counts contain a non-finite, negative or non-integer count')
        raw = sp.csr_matrix((values, cols, local_ptr), shape=(end - start, n_features), dtype=np.float64)
        raw.sum_duplicates()
        raw.eliminate_zeros()
        local_cell = source_to_cell[start:end]
        keep_local = np.flatnonzero(local_cell >= 0)
        if not len(keep_local):
            continue
        cell_rows = local_cell[keep_local]
        retained = raw[keep_local].tocsr()
        row_totals = np.asarray(retained.sum(axis=1)).ravel()
        recorded = full_totals[cell_rows]
        if (recorded <= 0).any() or not np.array_equal(row_totals, recorded):
            stop('raw full-source count row totals differ from pinned cell_metadata.full_count_total')
        selected_units = cell_unit[cell_rows]
        if score_matrices is None:
            active_wt = unit_active[selected_units] & (unit_genotype[selected_units] == 'WT')
            for population in ALLPOPS:
                choose = active_wt & (unit_pop[selected_units] == population)
                if not choose.any():
                    continue
                local_rows = np.flatnonzero(choose)
                block = retained[local_rows].tocsr()
                block.eliminate_zeros()
                factors = 10000.0 / recorded[local_rows]
                linear_sums[population] += np.asarray(block.multiply(factors[:, None]).sum(axis=0)).ravel()
                detections[population] += np.asarray(block.getnnz(axis=0)).ravel().astype(np.int64)
                wt_nuclei[population] += int(len(local_rows))
        if score_matrices is not None:
            is_target = unit_active[selected_units] & np.isin(unit_pop[selected_units], POPS)
            for population in POPS:
                choose = is_target & (unit_pop[selected_units] == population)
                if not choose.any():
                    continue
                local_rows = np.flatnonzero(choose)
                cell_ids = cell_rows[local_rows]
                unit_ids = selected_units[local_rows]
                local_ids = target_unit_maps[population][unit_ids]
                if (local_ids < 0).any():
                    stop('scored nucleus maps to a unit excluded from its per-mouse table')
                block = retained[local_rows].tocsr()
                block.eliminate_zeros()
                factors = 10000.0 / recorded[local_rows]
                row_id = np.repeat(np.arange(len(local_rows)), np.diff(block.indptr))
                log_block = sp.csr_matrix((np.log1p(block.data * factors[row_id]), block.indices.copy(), block.indptr.copy()), shape=block.shape)
                scored = log_block @ score_matrices[population]
                values = scored.toarray() if sp.issparse(scored) else np.asarray(scored)
                group = sp.csr_matrix(
                    (np.ones(len(local_ids), dtype=np.float64), (local_ids, np.arange(len(local_ids)))),
                    shape=(len(target_nuclei[population]), len(local_ids)),
                )
                target_sums[population] += np.asarray(group @ values)
                target_nuclei[population] += np.bincount(local_ids, minlength=len(target_nuclei[population]))


def compute_markers(registry, linear_sums, detections, wt_nuclei):
    n_features = len(registry)
    means = {name: linear_sums[name] / wt_nuclei[name] for name in ALLPOPS}
    prevalence = {name: detections[name] / wt_nuclei[name] for name in ALLPOPS}
    mark, candidates, ranks, marker_rows, set_rows = {}, {}, {}, [], []
    feature_index = np.arange(n_features, dtype=np.int64)
    for name in ALLPOPS:
        other_names = [n for n in ALLPOPS if n != name]
        max_other = np.max(np.vstack([means[n] for n in other_names]), axis=0)
        own = means[name]
        ratio = np.where(max_other > 0, own / np.maximum(max_other, np.finfo(float).tiny), np.where(own > 0, np.inf, np.nan))
        mean_ok = own >= MEAN_FLOOR_CP10K
        prev_ok = prevalence[name] >= PREVALENCE_FLOOR
        base = mean_ok & prev_ok
        mark[name], candidates[name], ranks[name] = {}, {}, {}
        for threshold in THRESH:
            eligible = np.flatnonzero(base & (ratio >= threshold))
            order = eligible[np.lexsort((eligible, -ratio[eligible]))]
            selected = order[:CAP]
            flag = np.zeros(n_features, dtype=bool)
            flag[selected] = True
            rank = np.zeros(n_features, dtype=np.int32)
            rank[order] = np.arange(1, len(order) + 1, dtype=np.int32)
            mark[name][threshold] = flag
            candidates[name][threshold] = base & (ratio >= threshold)
            ranks[name][threshold] = rank
            measurable = bool(len(selected))
            set_rows.append({
                'population': name, 'threshold': threshold,
                'threshold_role': 'PRIMARY' if threshold == PRIMARY else 'SECONDARY',
                'n_all_source_genes': n_features, 'n_WT_nuclei': int(wt_nuclei[name]),
                'n_candidates_after_expression_floor': int(candidates[name][threshold].sum()),
                'n_markers_after_top50_cap': int(len(selected)),
                'measurable': measurable,
                'measurability_status': 'measurable' if measurable else 'not_measurable_no_qualifying_markers',
                'expression_floor_CP10k': MEAN_FLOOR_CP10K,
                'detection_fraction_floor': PREVALENCE_FLOOR,
                'priority_measurability_target': name in MEASURABILITY_TARGETS,
            })
        frame = registry.copy()
        frame['population'] = name
        frame['n_WT_nuclei'] = int(wt_nuclei[name])
        frame['own_WT_mean_CP10k'] = own
        frame['max_other_WT_mean_CP10k'] = max_other
        frame['own_WT_nonzero_fraction'] = prevalence[name]
        frame['ratio_vs_max_other'] = ratio
        frame['passes_mean_floor'] = mean_ok
        frame['passes_detection_floor'] = prev_ok
        frame['candidate_after_floors'] = base
        for threshold in THRESH:
            key = tkey(threshold)
            frame['candidate_' + key] = candidates[name][threshold]
            frame['marker_' + key] = mark[name][threshold]
            frame['rank_' + key] = ranks[name][threshold]
        marker_rows.append(frame)
    return mark, candidates, ranks, pd.concat(marker_rows, ignore_index=True), pd.DataFrame(set_rows)


def build_score_matrices(mark):
    matrices, combos, meta = {}, {}, {}
    n_features = len(next(iter(mark.values()))[THRESH[0]])
    for scored in POPS:
        combo_list = [(scored, marker_pop, threshold) for marker_pop in ALLPOPS if marker_pop != scored for threshold in THRESH]
        rows, cols = [], []
        counts = {}
        for column, (_, marker_pop, threshold) in enumerate(combo_list):
            used = mark[marker_pop][threshold] & ~mark[scored][threshold]
            indices = np.flatnonzero(used)
            counts[(scored, marker_pop, threshold)] = int(len(indices))
            rows.extend(indices.tolist())
            cols.extend([column] * len(indices))
            for idx in indices:
                meta[(scored, marker_pop, threshold, int(idx))] = True
        matrices[scored] = sp.csc_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n_features, len(combo_list)))
        combos[scored] = combo_list
    return matrices, combos, meta


def run_analysis(base, h5ad_path, out, real):
    global CHECKS
    CHECKS = []
    started = time.time()
    cells, units, cell_unit, unit_pop, unit_genotype, unit_age, unit_n, unit_active = load_inputs(base, real)
    registry = load_gene_registry(h5ad_path, base + S6 + 'feature_registry.tsv')
    if len(registry) <= 1414:
        stop('the source marker universe must extend beyond the S06 1,414-feature panel')
    source_to_cell, source_shape = source_row_map(h5ad_path, cells, units, real)
    ok('full-source H5AD provides %d genes outside the S06 1,414-feature panel' % (len(registry) - 1414), len(registry))
    unit_rows = units.copy()
    unit_rows['genotype_norm'] = unit_genotype
    unit_rows['age_norm'] = unit_age
    unit_rows['cell_type_raw'] = unit_pop
    n_units = len(units)
    for name in ALLPOPS:
        observed = pd.read_csv(base + p7(name), sep='\t')
        expected_rows = np.flatnonzero((unit_pop[cell_unit] == name) & unit_active[cell_unit])
        expected = cells.iloc[expected_rows].reset_index(drop=True)
        if list(observed.columns) != list(cells.columns) or len(observed) != len(expected):
            stop('S07 nucleus_metadata row count/columns differ for ' + name)
        for col in ('obs_name', 'unit_index', 'full_count_total'):
            if not np.array_equal(observed[col].to_numpy(), expected[col].to_numpy()):
                stop('S06 cell rows do not equal pinned S07 nucleus_metadata for ' + name + ' (' + col + ')')
    ok('S06 retained all-cell row mapping reproduces pinned S07 nucleus_metadata for all 13 populations', len(ALLPOPS))
    linear_sums = {name: np.zeros(len(registry), dtype=np.float64) for name in ALLPOPS}
    detections = {name: np.zeros(len(registry), dtype=np.int64) for name in ALLPOPS}
    wt_nuclei = {name: 0 for name in ALLPOPS}
    with h5py.File(h5ad_path, 'r') as handle:
        iter_count_chunks(handle, source_shape, source_to_cell, cells, cell_unit, units, unit_pop, unit_genotype,
                          unit_n, unit_active, linear_sums, detections, wt_nuclei)
    for name in ALLPOPS:
        expected_n = int(((unit_pop[cell_unit] == name) & unit_active[cell_unit] & (unit_genotype[cell_unit] == 'WT')).sum())
        if wt_nuclei[name] != expected_n or expected_n == 0:
            stop('WT source row count differs from unit metadata for ' + name)
    ok('two-pass input preparation: WT pooled linear CP10k means and WT-nucleus detection for every full-source gene/population', len(ALLPOPS) * len(registry))
    mark, candidates, ranks, marker_genes, marker_sets = compute_markers(registry, linear_sums, detections, wt_nuclei)
    ok('pre-declared floors applied to all source genes: own mean >= %.6g CP10k and WT detection >= %.6g' % (MEAN_FLOOR_CP10K, PREVALENCE_FLOOR), len(ALLPOPS) * len(THRESH))
    ok('measurability report contains Granule and GABA-1 for every threshold, including explicit not-measurable statuses', 2 * len(THRESH))
    matrices, combos, _ = build_score_matrices(mark)
    target_unit_maps, target_sums, target_nuclei = {}, {}, {}
    for scored in POPS:
        ids = np.flatnonzero((unit_pop == scored) & unit_active)
        target_unit_maps[scored] = np.full(n_units, -1, dtype=np.int64)
        target_unit_maps[scored][ids] = np.arange(len(ids), dtype=np.int64)
        target_sums[scored] = np.zeros((len(ids), len(combos[scored])), dtype=np.float64)
        target_nuclei[scored] = np.zeros(len(ids), dtype=np.int64)
    with h5py.File(h5ad_path, 'r') as handle:
        iter_count_chunks(handle, source_shape, source_to_cell, cells, cell_unit, units, unit_pop, unit_genotype,
                          unit_n, unit_active, linear_sums, detections, wt_nuclei,
                          target_sums, target_nuclei, matrices, target_unit_maps)
    for scored in POPS:
        ids = np.flatnonzero((unit_pop == scored) & unit_active)
        if not np.array_equal(target_nuclei[scored], unit_n[ids]):
            stop('source nuclei do not reproduce active per-mouse totals for ' + scored)
    ok('all scored nuclei are assigned once to their per-mouse table; no pooling across ages', int(sum(len(x) for x in target_nuclei.values())))
    count_rows, mouse_rows, test_rows, summary_rows, marker_list_rows = [], [], [], [], []
    for scored in POPS:
        ids = np.flatnonzero((unit_pop == scored) & unit_active)
        ordered = units.iloc[ids].copy()
        ordered['g'] = unit_genotype[ids]
        ordered['age'] = unit_age[ids]
        ordered = ordered.sort_values(['age', 'g', 'animal_id'])
        local_by_unit = target_unit_maps[scored]
        for col, (_, marker_pop, threshold) in enumerate(combos[scored]):
            marker_before = mark[marker_pop][threshold]
            used = marker_before & ~mark[scored][threshold]
            n_used = int(used.sum())
            n_candidates = int(candidates[marker_pop][threshold].sum())
            n_selected = int(marker_before.sum())
            n_self = int((marker_before & mark[scored][threshold]).sum())
            count_rows.append({
                'scored_population': scored, 'marker_population': marker_pop, 'threshold': threshold,
                'threshold_role': 'PRIMARY' if threshold == PRIMARY else 'SECONDARY',
                'n_candidates_pre_cap': n_candidates, 'n_markers_after_cap': n_selected,
                'n_capped_by_ratio': n_candidates - n_selected,
                'n_removed_markers_of_scored_population': n_self, 'n_markers_used': n_used,
                'status': 'scored' if n_used else ('no_markers_after_self_removal' if n_selected else 'no_qualifying_markers'),
            })
            for feature_index in np.flatnonzero(used):
                marker_list_rows.append({
                    'scored_population': scored, 'marker_population': marker_pop, 'threshold': threshold,
                    'feature_index': int(feature_index), 'stable_id': registry['stable_id'].iloc[feature_index],
                    'dataset_symbol': registry['dataset_symbol'].iloc[feature_index],
                    'own_WT_mean_CP10k': float(linear_sums[marker_pop][feature_index] / wt_nuclei[marker_pop]),
                    'own_WT_nonzero_fraction': float(detections[marker_pop][feature_index] / wt_nuclei[marker_pop]),
                })
            values = target_sums[scored][:, col] / target_nuclei[scored]
            if n_used:
                for local_id in range(len(ids)):
                    row = ordered.iloc[np.flatnonzero(local_by_unit[ordered['unit_index'].to_numpy()] == local_id)[0]]
                    mouse_rows.append({
                        'scored_population': scored, 'marker_population': marker_pop, 'threshold': threshold,
                        'threshold_role': 'PRIMARY' if threshold == PRIMARY else 'SECONDARY',
                        'animal_id': str(row['animal_id']), 'genotype': row['g'], 'age_weeks': int(row['age']),
                        'n_nuclei': int(unit_n[int(row['unit_index'])]), 'n_markers_used': n_used,
                        'mean_score': float(values[local_id]),
                    })
            for age in AGES:
                wt_units = ordered[(ordered['age'] == age) & (ordered['g'] == 'WT')]
                ki_units = ordered[(ordered['age'] == age) & (ordered['g'] == 'KI')]
                wt_ids = local_by_unit[wt_units['unit_index'].to_numpy(dtype=np.int64)]
                ki_ids = local_by_unit[ki_units['unit_index'].to_numpy(dtype=np.int64)]
                wt_values = values[wt_ids] if n_used else np.asarray([], dtype=float)
                ki_values = values[ki_ids] if n_used else np.asarray([], dtype=float)
                record = {
                    'scored_population': scored, 'marker_population': marker_pop, 'threshold': threshold,
                    'threshold_role': 'PRIMARY' if threshold == PRIMARY else 'SECONDARY', 'age_weeks': age,
                    'n_markers_used': n_used, 'n_WT': int(len(wt_units)), 'n_KI': int(len(ki_units)),
                    'status': 'scored' if n_used else ('no_markers_after_self_removal' if n_selected else 'no_qualifying_markers'),
                }
                if n_used and len(wt_values) == 3 and len(ki_values) == 3:
                    wt_mean, ki_mean = float(wt_values.mean()), float(ki_values.mean())
                    diff = ki_mean - wt_mean
                    pct = 100.0 * diff / wt_mean if wt_mean != 0 else np.nan
                    record.update({
                        'WT_mean_score': wt_mean, 'KI_mean_score': ki_mean,
                        'diff_KI_minus_WT': diff, 'diff_pct_of_WT_mean_score': pct,
                        'pct_status': 'defined' if wt_mean != 0 else 'undefined_WT_mean_zero',
                        'n_labelings': 20,
                        'WT_values': ';'.join('%s:%.10g' % (str(row['animal_id']), float(values[local_by_unit[int(row['unit_index'])]])) for _, row in wt_units.iterrows()),
                        'KI_values': ';'.join('%s:%.10g' % (str(row['animal_id']), float(values[local_by_unit[int(row['unit_index'])]])) for _, row in ki_units.iterrows()),
                    })
                    record.update(exact_perm(wt_values, ki_values))
                elif n_used:
                    record['status'] = 'not_evaluable_not_3v3'
                test_rows.append(record)
                summary_rows.append({
                    'scored_population': scored, 'marker_population': marker_pop, 'threshold': threshold,
                    'threshold_role': 'PRIMARY' if threshold == PRIMARY else 'SECONDARY', 'age_weeks': age,
                    'n_markers_used': n_used, 'status': record['status'],
                    'WT_mean_score': record.get('WT_mean_score', np.nan),
                    'diff_KI_minus_WT': record.get('diff_KI_minus_WT', np.nan),
                    'diff_pct_of_WT_mean_score': record.get('diff_pct_of_WT_mean_score', np.nan),
                    'pct_status': record.get('pct_status', 'not_scored'),
                    'p_two_sided': record.get('p_two_sided', np.nan),
                })
    tests = pd.DataFrame(test_rows)
    if len(tests) != len(POPS) * (len(ALLPOPS) - 1) * len(THRESH) * len(AGES):
        stop('within-age test table does not contain every scored-population/marker-population/threshold/age combination')
    scored_tests = tests[tests['n_markers_used'] > 0]
    if len(scored_tests) and not (scored_tests['n_labelings'].fillna(0).astype(int) == 20).all():
        stop('a scored per-age KI-vs-WT comparison lacks all 20 exact labelings')
    ok('per-mouse KI-vs-WT table has all %d within-age contrasts; scored contrasts use 3 WT, 3 KI and 20 exact labelings' % len(tests), len(tests))
    measurability = marker_sets.copy()
    measurability['priority_status'] = np.where(
        measurability['priority_measurability_target'] & measurability['measurable'], 'priority_target_measurable',
        np.where(measurability['priority_measurability_target'], 'PRIORITY_TARGET_NOT_MEASURABLE', 'not_priority_target'))
    os.makedirs(out, exist_ok=True)
    if os.listdir(out):
        stop('the output directory must exist and be empty')

    def write_table(frame, name):
        frame.to_csv(out + name, sep='\t', index=False, float_format='%.10g', na_rep='NA', mode='x')

    write_table(marker_genes, 'marker_genes.tsv')
    write_table(measurability, 'marker_set_measurability.tsv')
    write_table(pd.DataFrame(count_rows), 'marker_counts.tsv')
    write_table(pd.DataFrame(marker_list_rows, columns=['scored_population', 'marker_population', 'threshold', 'feature_index', 'stable_id', 'dataset_symbol', 'own_WT_mean_CP10k', 'own_WT_nonzero_fraction']), 'marker_lists.tsv')
    write_table(pd.DataFrame(mouse_rows, columns=['scored_population', 'marker_population', 'threshold', 'threshold_role', 'animal_id', 'genotype', 'age_weeks', 'n_nuclei', 'n_markers_used', 'mean_score']), 'marker_mouse_means.tsv')
    write_table(tests, 'marker_kivswt_within_age.tsv')
    write_table(pd.DataFrame(summary_rows), 'marker_summary.tsv')
    write_table(pd.DataFrame(CHECKS, columns=['check_passed', 'items']), 'checks.tsv')
    unmeasurable = measurability[~measurability['measurable']]
    priority_unmeasurable = unmeasurable[unmeasurable['priority_measurability_target']]
    result = {
        'mode': 'real' if real else 'fixture', 'checks': len(CHECKS),
        'source_h5ad_sha256': H5AD_SHA256 if real else sha256(h5ad_path),
        'source_dimensions': {'n_obs': int(source_shape[0]), 'n_all_source_genes': int(source_shape[1]), 's06_panel_genes': 1414},
        'populations': ALLPOPS, 'scored_populations': POPS, 'ages_weeks': AGES,
        'thresholds': THRESH, 'primary_threshold': PRIMARY,
        'marker_universe': 'every gene column in the pinned source H5AD layers/counts matrix and var gene_ids; no 1,414-feature-panel restriction',
        'expression_rule_declared_before_run': {
            'scale': RATIO_SCALE, 'own_WT_mean_CP10k_minimum_inclusive': MEAN_FLOOR_CP10K,
            'own_WT_nucleus_detection_fraction_minimum_inclusive': PREVALENCE_FLOOR,
            'ratio_against_each_other_population_minimum_inclusive': THRESH,
            'cap': CAP, 'cap_order': 'descending pooled WT ratio; ascending full-source feature_index tie-break',
        },
        'measurability_targets': MEASURABILITY_TARGETS,
        'unmeasurable_marker_sets': unmeasurable[['population', 'threshold', 'measurability_status']].to_dict(orient='records'),
        'priority_target_unmeasurable_marker_sets': priority_unmeasurable[['population', 'threshold', 'measurability_status']].to_dict(orient='records'),
        'per_mouse_comparisons': len(tests), 'scored_within_age_tests': int(len(scored_tests)),
        'exact_labelings_per_scored_test': 20,
        'score_scale': SCORE_SCALE,
        'description': 'per mouse, KI vs WT within each age; exact 3-v-3 label permutations; descriptive, no multiplicity adjustment',
        'caveat': 'a KI-vs-WT difference in a neighbour-marker score can be leak-in OR real expression; T6b cannot distinguish the two',
        'checks_summary': CHECKS, 'python': sys.version.split()[0], 'numpy': np.__version__,
        'pandas': pd.__version__, 'scipy': scipy.__version__,
        'analysis_script_sha256': sha256(os.path.abspath(__file__)), 'seconds': round(time.time() - started, 2),
    }
    with open(out + 'RESULT.json', 'x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(out)


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--real', action='store_true')
    group.add_argument('--fixture')
    parser.add_argument('out')
    args = parser.parse_args()
    if not os.path.isdir(args.out) or os.listdir(args.out):
        stop('the output directory must exist and be empty')
    if args.real:
        base, pins = REAL_CFG['root'], REAL_CFG['pins']
        h5ad_path = base + H5AD_REL
    else:
        base = args.fixture.rstrip('/') + '/'
        pins = json.load(open(base + 'FIXTURE_PINS.json', encoding='utf-8'))
        h5ad_path = base + H5AD_REL
    verify_inputs(base, pins)
    run_analysis(base, h5ad_path, args.out.rstrip('/') + '/', args.real)


if __name__ == '__main__':
    main()
