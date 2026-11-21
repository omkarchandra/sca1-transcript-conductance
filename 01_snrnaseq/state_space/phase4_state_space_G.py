"""All-gene G(a)/G(b) builder for the shared P/G state-space interface.

Draft-only real-data entry point. It opens the pinned H5AD only after a
G_builder-scoped, byte-bound review review and this run's <AUTHOR> yes pass.
Output matrices are float64
log1p(CP10k), with aligned nucleus and feature tables.
"""
from __future__ import annotations

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)


import argparse
import hashlib
import json
import os
import sys
import time
import traceback
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
try:
    from anndata.io import read_elem
except ModuleNotFoundError:
    # AnnData 0.10.x exposes this reader from the private namespace; 0.11+ also
    # provides the public anndata.io import. cluster's single_cell env is 0.10.8.
    from anndata._io.specs import read_elem

HERE = Path(__file__).resolve().parent
ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
sys.path.insert(0, str(HERE))
import phase4_state_space_shared_A as shared

POPS = ('Purkinje cell', 'Bergmann glia')
T5_DISCOVERY_POPULATIONS = (
    'Astrocyte', 'Bergmann glia', 'Endothelial cell', 'GABAergic interneuron 1',
    'GABAergic interneuron 2', 'GABAergic interneuron 3', 'Granule cell', 'Microglia',
    'Oligodendrocyte', 'Oligodendrocyte progenitor cell', 'Pericyte', 'Purkinje cell',
    'Unipolar brush cell',
)
POP_SLUG = {p: p.replace(' ', '_') for p in POPS}
DISCOVERY_AGES = (5, 18, 24, 30)
READOUT_AGES = (5, 12, 18, 24, 30)
TWELVE_MICE = ('2061', '2062', '22018', '2063', '2065', '22019')
H5AD_SHA256 = '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'
NUCLEUS_MANIFEST_SHA256 = 'c3e7d228e97c4991b7d288f46e210308b3ace0cb839aac386f3d578730d08352'
ANIMAL_MANIFEST_SHA256 = '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'
SCOPING_RECORD_SHA256 = 'a7131e744951ec404f15ba003542a80ad793df814ea17c58f6290a05b8ec7447'
TWELVE_FULL_RECORD_SHA256 = '6b0e3c442b1e6ed41ecd2d67ee1d49e59149c22cfcafedadaab7956f0b784123'
AMBIENT_RELEASE_SHA256 = '40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298'
AMBIENT_RELEASE_REL = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz')
T5_RATIO_SCRIPT_REL = _release_paths.expand('${REPO_ROOT}/01_snrnaseq/s12_own_gene/s12_t5_decomp_r1.py')
T5_RATIO_SCRIPT_SHA256 = 'a2484ba3afe60830644d4751390db9bda66d6771472c6f0640a45ac8d23db29f'
SCOPING_REL = _release_paths.expand('${DATA_ROOT}/phase4_12w_scoping/SCOPING.json')
TWELVE_FULL_RECORD_REL = _release_paths.expand('${DATA_ROOT}/phase4_12w_full/RECORD.json')
CHUNK_ROWS = 1024
WT_UNIT_FLOOR = 20

# All records are manifest rows for either target population. Manifest included=False
# for the six 12 W mice has exactly one reason (the design's excluded age); no second
# nucleus-level filter is imposed. The two cited Phase-4 records lock this re-inclusion.
TWELVE_QC_POLICY = {
    'record_scoping_sha256': SCOPING_RECORD_SHA256,
    'record_full_process_sha256': TWELVE_FULL_RECORD_SHA256,
    'manifest_rule': 're-include 12 W manifest rows only when their sole exclusion_reason is age_weeks in design.excluded_ages',
    'population_rule': 'retain all 12 W rows whose cell_type_raw is Purkinje cell or Bergmann glia; do not apply is_target because it is Purkinje-specific in this manifest',
    'additional_nucleus_filters': [],
}


def sha(path: Path | str) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _strict_boolean_labels(series: pd.Series, column: str) -> np.ndarray:
    labels = series.astype(str).str.strip().str.upper()
    unexpected = sorted(set(labels) - {'TRUE', 'FALSE'})
    assert not unexpected, (column, 'expected only TRUE/FALSE labels', unexpected)
    return labels.eq('TRUE').to_numpy()


def t5_pooled_ratio_of(own: np.ndarray, max_other: np.ndarray) -> np.ndarray:
    """Port of S12 T5 main().ratio_of (pinned script, lines 214-215)."""
    own = np.asarray(own, dtype=np.float64)
    max_other = np.asarray(max_other, dtype=np.float64)
    assert own.shape == max_other.shape and np.isfinite(own).all() and np.isfinite(max_other).all()
    assert (own >= 0).all() and (max_other >= 0).all()
    ratio = np.divide(own, max_other, out=np.zeros_like(own), where=max_other > 0)
    ratio[(max_other == 0) & (own > 0)] = np.inf
    ratio[(max_other == 0) & (own == 0)] = np.nan
    return ratio


def check_input_pins(h5ad: Path, manifest: Path, animal_manifest: Path,
                     scoping_record: Path, full_record: Path, ambient_release: Path,
                     t5_ratio_script: Path) -> dict:
    expected = {
        'h5ad': (h5ad, H5AD_SHA256),
        'nucleus_manifest': (manifest, NUCLEUS_MANIFEST_SHA256),
        'animal_manifest': (animal_manifest, ANIMAL_MANIFEST_SHA256),
        'scoping_record': (scoping_record, SCOPING_RECORD_SHA256),
        'phase4_12w_full_record': (full_record, TWELVE_FULL_RECORD_SHA256),
        'S08_S11_ambient_release': (ambient_release, AMBIENT_RELEASE_SHA256),
        'S12_T5_ratio_reference_script': (t5_ratio_script, T5_RATIO_SCRIPT_SHA256),
    }
    got = {}
    for label, (path, digest) in expected.items():
        actual = sha(path)
        assert actual == digest, (label, str(path), actual, digest)
        got[str(path)] = actual
    return got


def _manifest_bool(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.strip().str.lower().eq('true').to_numpy()


def build_manifest_selection(manifest: pd.DataFrame, animals: pd.DataFrame,
                             scoping: dict, strict_expected_counts: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply frozen discovery inclusion plus explicit 12 W age-only re-inclusion."""
    required = {'obs_name', 'animal_id', 'cell_type_raw', 'included', 'exclusion_reason',
                'age_weeks_obs', 'genotype_final'}
    assert required <= set(manifest.columns), required - set(manifest.columns)
    assert {'animal_id', 'age_weeks', 'genotype', 'sex'} <= set(animals.columns)
    assert not manifest.obs_name.duplicated().any()
    assert not animals.animal_id.astype(str).duplicated().any()
    manifest = manifest.copy()
    manifest['animal_id'] = manifest.animal_id.astype(str)
    animals = animals.copy()
    animals['animal_id'] = animals.animal_id.astype(str)
    animal_lookup = animals.set_index('animal_id')
    assert set(TWELVE_MICE) <= set(animal_lookup.index)
    expected_six = scoping['six_mice']
    assert set(expected_six) == set(TWELVE_MICE)
    # These six mice were excluded only by age in the locked manifest. Fail closed if
    # an unrelated QC exclusion appears rather than silently admitting or dropping it.
    twelve = manifest[manifest.animal_id.isin(TWELVE_MICE)]
    if strict_expected_counts:
        assert len(twelve) == 38017
    assert (twelve.included.astype(str).str.lower() == 'false').all()
    assert twelve.exclusion_reason.eq('age_weeks in design.excluded_ages').all()
    assert twelve.age_weeks_obs.astype(int).eq(12).all()
    assert (twelve.genotype_final.astype(str).str.upper().isin({'WT', 'KI'})).all()
    for animal_id, row in expected_six.items():
        assert int((twelve.animal_id == animal_id).sum()) == int(row['n_nuclei_table_s1'])
        am = animal_lookup.loc[animal_id]
        assert str(am.genotype).upper() == str(row['genotype']).upper()
        assert int(am.age_weeks) == 12
        assert str(am.sex).upper() == str(row['sex']).upper()

    included = _manifest_bool(manifest.included)
    discovery = included & manifest.age_weeks_obs.astype(int).isin(DISCOVERY_AGES).to_numpy()
    twelve_keep = manifest.animal_id.isin(TWELVE_MICE).to_numpy()
    target = manifest.cell_type_raw.isin(POPS).to_numpy()
    # All included discovery populations are needed for the G(b) own/max-other ratio;
    # only the two target populations' 12 W rows are admitted to output matrices.
    selected = discovery | (twelve_keep & target)
    assert selected.any()
    assert manifest.loc[selected, 'age_weeks_obs'].astype(int).isin(READOUT_AGES).all()
    selected_rows = manifest.loc[selected].copy().reset_index(drop=True)
    selected_rows = selected_rows.merge(
        animals[['animal_id', 'age_weeks', 'genotype', 'sex']], on='animal_id', how='left',
        validate='many_to_one', suffixes=('', '_animal'))
    assert selected_rows.age_weeks.notna().all() and selected_rows.genotype.notna().all() and selected_rows.sex.notna().all()
    assert (selected_rows.age_weeks.astype(int) == selected_rows.age_weeks_obs.astype(int)).all()
    assert (selected_rows.genotype.astype(str).str.upper() == selected_rows.genotype_final.astype(str).str.upper()).all()
    selected_rows['age_weeks'] = selected_rows.age_weeks.astype(int)
    selected_rows['genotype'] = selected_rows.genotype.astype(str).str.upper()
    selected_rows['sex'] = selected_rows.sex.astype(str).str.upper()
    assert selected_rows.genotype.isin(['WT', 'KI']).all()
    assert not selected_rows.obs_name.duplicated().any()

    # The per-mouse retained count is compared against the accepted scoping record.
    retained = selected_rows[selected_rows.animal_id.isin(TWELVE_MICE)].groupby(
        ['animal_id', 'cell_type_raw'], sort=True).size()
    count_rows = []
    for animal_id in TWELVE_MICE:
        expected_types = scoping['per_mouse_by_celltype_counts'][animal_id]
        for population in POPS:
            expected_n = int(expected_types[population]['nuclei'])
            kept_n = int(retained.get((animal_id, population), 0))
            assert kept_n == expected_n, (animal_id, population, kept_n, expected_n)
            count_rows.append(dict(mouse_id=animal_id,
                                   genotype=str(animal_lookup.loc[animal_id, 'genotype']).upper(),
                                   age_weeks=12, population=population,
                                   manifest_nuclei=expected_n, retained_nuclei=kept_n,
                                   manifest_included_true=int(expected_types[population]['included']),
                                   reinclude_reason='age-only; §14 scoping contract', silently_dropped=0))
    if strict_expected_counts:
        assert sum(row['retained_nuclei'] for row in count_rows) == 1144
    return selected_rows, pd.DataFrame(count_rows)


def ambient_flagged_ids(path: Path) -> dict[str, set[str]]:
    table = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False,
                        usecols=['stable_id', 'population_registry_name',
                                 'potential_ambient_suspect', 'ambient_flagged_hit'])
    result = {}
    for population in POPS:
        rows = table[table.population_registry_name == population]
        flagged = (_strict_boolean_labels(rows.potential_ambient_suspect, 'potential_ambient_suspect')
                   | _strict_boolean_labels(rows.ambient_flagged_hit, 'ambient_flagged_hit'))
        result[population] = set(rows.loc[flagged, 'stable_id'].astype(str))
        assert len(rows) and len(result[population])
    return result


def _iter_chunks(n_rows: int):
    for start in range(0, n_rows, CHUNK_ROWS):
        yield start, min(start + CHUNK_ROWS, n_rows)


def _build_outputs(h5ad: Path, selected: pd.DataFrame, var: pd.DataFrame,
                   ambient_ids: dict[str, set[str]], out_dir: Path,
                   qc_counts: pd.DataFrame, expected_manifest_obs=None) -> dict:
    """Two-pass CSR read: choose features, then write aligned dense float64 matrices."""
    if out_dir.exists():
        raise FileExistsError(f'refusing to overwrite existing output: {out_dir}')
    out_dir.mkdir(parents=True)
    names = var.index.astype(str).to_numpy()
    stable_ids = var['gene_ids'].astype(str).to_numpy()
    assert len(names) == len(stable_ids) and not pd.Index(stable_ids).duplicated().any()
    n_genes = len(var)
    selected_lookup = selected.set_index('obs_name', verify_integrity=True)
    WT_DISCOVERY = selected.genotype.eq('WT') & selected.age_weeks.isin(DISCOVERY_AGES)
    selected_rows_by_pop = {
        pop: selected.loc[selected.cell_type_raw.eq(pop)].copy().reset_index(drop=True)
        for pop in POPS
    }
    for pop, frame in selected_rows_by_pop.items():
        assert set(frame.age_weeks.astype(int)) == set(READOUT_AGES), (pop, sorted(frame.age_weeks.unique()))
        assert set(frame.genotype) == {'WT', 'KI'}
        assert not frame.obs_name.duplicated().any()
    target_names = set(POPS)
    selected_obs = set(selected.obs_name.astype(str))

    # Match S12 T5's pooled-WT denominator: mouse x cell-type units with at least
    # 20 WT discovery nuclei; this is a mean-coverage rule, never an output filter.
    unit_sizes = selected.loc[WT_DISCOVERY].groupby(['animal_id', 'cell_type_raw']).size()
    wt_unit_ok = {(str(a), str(p)) for (a, p), n in unit_sizes.items() if int(n) >= WT_UNIT_FLOOR}
    all_types = tuple(sorted(selected.loc[WT_DISCOVERY, 'cell_type_raw'].unique()))
    if expected_manifest_obs is not None:
        assert set(all_types) == set(T5_DISCOVERY_POPULATIONS), (
            'T5 denominator population universe changed', sorted(all_types))
    sums_by_type = {p: np.zeros(n_genes, dtype=np.float64) for p in all_types}
    n_by_type = {p: 0 for p in all_types}
    detected_by_target = {p: np.zeros(n_genes, dtype=np.int64) for p in POPS}
    n_detected_by_target = {p: 0 for p in POPS}
    target_matrix_rows = {p: np.zeros(len(selected_rows_by_pop[p]), dtype=bool) for p in POPS}
    target_local_indices = {p: {} for p in POPS}
    for pop, frame in selected_rows_by_pop.items():
        target_local_indices[pop] = {str(name): int(i) for i, name in enumerate(frame.obs_name)}

    with h5py.File(h5ad, 'r') as f:
        obs = read_elem(f['obs'])
        assert obs.index.is_unique and len(obs) == f['layers/counts/indptr'].shape[0] - 1
        obs_names = obs.index.astype(str).to_numpy()
        if expected_manifest_obs is not None:
            assert len(obs_names) == len(expected_manifest_obs) and set(obs_names) == expected_manifest_obs
        assert selected_obs <= set(obs_names)
        selected_by_obs = selected_lookup.reindex(obs_names)
        selected_by_obs['obs_name'] = obs_names
        selected_mask = selected_by_obs.animal_id.notna().to_numpy()
        assert int(selected_mask.sum()) == len(selected)
        counts_group = f['layers/counts']
        indptr = counts_group['indptr'][...].astype(np.int64)
        for start, stop in _iter_chunks(len(obs_names)):
            lo, hi = int(indptr[start]), int(indptr[stop])
            values = counts_group['data'][lo:hi].astype(np.float64)
            columns = counts_group['indices'][lo:hi].astype(np.int64)
            row_idx = np.repeat(np.arange(stop - start), np.diff(indptr[start:stop + 1]))
            raw = sp.csr_matrix((values, (row_idx, columns)), shape=(stop - start, n_genes))
            raw.sum_duplicates(); raw.eliminate_zeros()
            totals = np.asarray(raw.sum(axis=1)).ravel()
            assert np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all()
            global_indices = np.arange(start, stop, dtype=np.int64)
            selected_here = selected_mask[start:stop]
            if not selected_here.any():
                continue
            local_sel = np.flatnonzero(selected_here)
            sel_meta = selected_by_obs.iloc[start + local_sel]
            full_rows = raw[local_sel]
            sel_totals = totals[local_sel]
            assert (sel_totals > 0).all(), 'selected nucleus has zero full-transcriptome UMI'
            cp10k = full_rows.multiply((10000.0 / sel_totals)[:, None]).tocsr()
            discovery_wt = (sel_meta.genotype.astype(str).str.upper().eq('WT').to_numpy()
                            & sel_meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy())
            if discovery_wt.any():
                wt_meta = sel_meta.iloc[np.flatnonzero(discovery_wt)]
                wt_cp = cp10k[np.flatnonzero(discovery_wt)]
                for cell_type in all_types:
                    type_rows = np.flatnonzero(wt_meta.cell_type_raw.to_numpy() == cell_type)
                    if not len(type_rows):
                        continue
                    if cell_type in target_names:
                        own_sum = np.asarray(wt_cp[type_rows].sum(axis=0)).ravel()
                        detected_by_target[cell_type] += np.asarray((full_rows[np.flatnonzero(discovery_wt)[type_rows]] > 0).sum(axis=0)).ravel().astype(np.int64)
                        n_detected_by_target[cell_type] += int(len(type_rows))
                    animal_values = wt_meta.iloc[type_rows].animal_id.astype(str).to_numpy()
                    unit_rows = np.array([(animal, cell_type) in wt_unit_ok for animal in animal_values], dtype=bool)
                    if unit_rows.any():
                        sums_by_type[cell_type] += np.asarray(wt_cp[type_rows[unit_rows]].sum(axis=0)).ravel()
                        n_by_type[cell_type] += int(unit_rows.sum())
            # Store source row positions in output row order; this also proves every eligible row was visited.
            for pop in POPS:
                pop_mask = sel_meta.cell_type_raw.astype(str).eq(pop).to_numpy()
                if pop_mask.any():
                    obs_for_pop = sel_meta.iloc[np.flatnonzero(pop_mask)].obs_name.astype(str).to_numpy()
                    for obs_name in obs_for_pop:
                        target_matrix_rows[pop][target_local_indices[pop][obs_name]] = True
        for pop in POPS:
            missing = np.flatnonzero(~target_matrix_rows[pop])
            assert not len(missing), (f'{pop}: eligible nuclei were not visited in H5AD pass 1',
                                      selected_rows_by_pop[pop].iloc[missing[:5]].obs_name.tolist(), len(missing))

    assert all(n_by_type[p] > 0 for p in all_types), 'a WT population has no >=20-nucleus units'
    means_by_type = {
        p: (sums_by_type[p] / n_by_type[p] if n_by_type[p] else np.zeros(n_genes, dtype=np.float64))
        for p in all_types
    }
    features_by_pop: dict[str, dict[str, np.ndarray]] = {}
    feature_metrics: dict[str, pd.DataFrame] = {}
    for pop in POPS:
        wt_rows = selected.loc[WT_DISCOVERY & selected.cell_type_raw.eq(pop)]
        assert len(wt_rows) and n_detected_by_target[pop] == len(wt_rows)
        assert n_by_type.get(pop, 0) > 0, f'{pop}: no WT discovery unit meets the 20-nucleus mean floor'
        gene_detection_fraction = detected_by_target[pop] / n_detected_by_target[pop]
        ga = np.flatnonzero(gene_detection_fraction >= .01)
        assert len(ga), f'{pop}: G(a) gene set is empty'
        other = [p for p in all_types if p != pop]
        max_other = np.max(np.vstack([means_by_type[p] for p in other]), axis=0)
        own = means_by_type[pop]
        # Same own/max-other edge semantics as S12 T5's main().ratio_of, now
        # applied to linear CP10k means over the all-gene universe.
        ratio = t5_pooled_ratio_of(own, max_other)
        own_detected_fraction = gene_detection_fraction
        gb_mask = (ratio >= 1.0) & (own >= 1.0) & (own_detected_fraction >= .20)
        flagged = np.array([stable_id in ambient_ids[pop] for stable_id in stable_ids], dtype=bool)
        gb_mask &= ~flagged
        gb = np.flatnonzero(gb_mask)
        assert len(gb), f'{pop}: G(b) gene set is empty'
        features_by_pop[pop] = {'G_a': ga, 'G_b': gb}
        feature_metrics[pop] = pd.DataFrame({
            'feature_index': np.arange(n_genes, dtype=int), 'stable_id': stable_ids,
            'symbol': names, 'WT_mean_CP10k_own': own, 'WT_max_other_CP10k': max_other,
            'WT_ratio_own_over_max_other': ratio, 'WT_detection_fraction_own': own_detected_fraction,
            'S08_ambient_flagged': flagged,
            'gene_selection_population': 'WT discovery nuclei only; ages 5,18,24,30; KI and age 12 excluded',
            'selected_G_a': gene_detection_fraction >= .01,
            'selected_G_b': gb_mask,
        })

    matrices = {}
    metadata_frames = {}
    metadata_values = {}
    for pop in POPS:
        frame = selected_rows_by_pop[pop].reset_index(drop=True)
        metadata_frames[pop] = frame
        metadata_values[pop] = dict(log10_total_umi=np.full(len(frame), np.nan, dtype=np.float64),
                                    n_detected=np.full(len(frame), -1, dtype=np.int64))
        for space in ('G_a', 'G_b'):
            feature_idx = features_by_pop[pop][space]
            tag = f'{space}_{POP_SLUG[pop]}'
            matrix_path = out_dir / f'state_{tag}.npy'
            matrices[(pop, space)] = np.lib.format.open_memmap(
                matrix_path, mode='w+', dtype=np.float64, shape=(len(frame), len(feature_idx)))
            features = pd.DataFrame({'feature_index': np.arange(len(feature_idx), dtype=int),
                                     'stable_id': stable_ids[feature_idx], 'symbol': names[feature_idx],
                                     'var_position': feature_idx})
            features.to_csv(out_dir / f'state_{tag}_features.tsv', sep='\t', index=False)
    # Pass two writes each row once to both requested gene spaces using the complete-transcriptome total.
    with h5py.File(h5ad, 'r') as f:
        obs = read_elem(f['obs'])
        obs_names = obs.index.astype(str).to_numpy()
        selected_by_obs = selected_lookup.reindex(obs_names)
        selected_by_obs['obs_name'] = obs_names
        selected_mask = selected_by_obs.animal_id.notna().to_numpy()
        counts_group = f['layers/counts']
        indptr = counts_group['indptr'][...].astype(np.int64)
        for start, stop in _iter_chunks(len(obs_names)):
            lo, hi = int(indptr[start]), int(indptr[stop])
            values = counts_group['data'][lo:hi].astype(np.float64)
            columns = counts_group['indices'][lo:hi].astype(np.int64)
            row_idx = np.repeat(np.arange(stop - start), np.diff(indptr[start:stop + 1]))
            raw = sp.csr_matrix((values, (row_idx, columns)), shape=(stop - start, n_genes))
            raw.sum_duplicates(); raw.eliminate_zeros()
            totals = np.asarray(raw.sum(axis=1)).ravel()
            selected_here = selected_mask[start:stop]
            if not selected_here.any():
                continue
            local_rows = np.flatnonzero(selected_here)
            sel_meta = selected_by_obs.iloc[start + local_rows]
            selected_counts = raw[local_rows]
            selected_totals = totals[local_rows]
            cp10k = selected_counts.multiply((10000.0 / selected_totals)[:, None]).tocsr()
            cp10k.data = np.log1p(cp10k.data)
            for pop in POPS:
                pop_positions = np.flatnonzero(sel_meta.cell_type_raw.astype(str).eq(pop).to_numpy())
                if not len(pop_positions):
                    continue
                obs_for_pop = sel_meta.iloc[pop_positions].obs_name.astype(str).to_numpy()
                out_rows = np.array([target_local_indices[pop][name] for name in obs_for_pop], dtype=int)
                totals_pop = selected_totals[pop_positions]
                full_counts_pop = selected_counts[pop_positions]
                n_detected = np.asarray((full_counts_pop > 0).sum(axis=1)).ravel().astype(np.int64)
                metadata_values[pop]['log10_total_umi'][out_rows] = np.log10(totals_pop)
                metadata_values[pop]['n_detected'][out_rows] = n_detected
                for space in ('G_a', 'G_b'):
                    feature_idx = features_by_pop[pop][space]
                    selected_matrix = cp10k[pop_positions][:, feature_idx].toarray().astype(np.float64, copy=False)
                    matrices[(pop, space)][out_rows] = selected_matrix
        for pop in POPS:
            frame = metadata_frames[pop]
            log_total = metadata_values[pop]['log10_total_umi']
            n_detected = metadata_values[pop]['n_detected']
            assert np.isfinite(log_total).all() and (n_detected >= 0).all()
            metadata_out = pd.DataFrame({
                'nucleus_id': frame.obs_name.astype(str).to_numpy(),
                'mouse_id': frame.animal_id.astype(str).to_numpy(),
                'age_weeks': frame.age_weeks.to_numpy(dtype=int),
                'genotype': frame.genotype.astype(str).to_numpy(),
                'sex': frame.sex.astype(str).to_numpy(),
                'log10_total_umi': log_total,
                'n_detected': n_detected,
            })
            for space in ('G_a', 'G_b'):
                tag = f'{space}_{POP_SLUG[pop]}'
                metadata_out.to_csv(out_dir / f'state_{tag}_nuclei.tsv', sep='\t', index=False)
    for mmap in matrices.values():
        mmap.flush()
    del matrices
    qc_counts.to_csv(out_dir / 'TWELVE_WEEK_RETAINED_COUNTS.tsv', sep='\t', index=False)
    feature_dir = out_dir / 'FEATURE_SELECTION'
    feature_dir.mkdir()
    for pop in POPS:
        feature_metrics[pop].to_csv(feature_dir / f'{POP_SLUG[pop]}_GENE_SET_METRICS.tsv.gz',
                                    sep='\t', index=False, compression='gzip')
    return {
        'outputs': sorted(p.name for p in out_dir.iterdir() if p.is_file()),
        'n_output_nuclei': {pop: int(len(metadata_frames[pop])) for pop in POPS},
        'n_features': {pop: {space: int(len(features_by_pop[pop][space])) for space in ('G_a', 'G_b')} for pop in POPS},
        '12w_retained_total': int(qc_counts.retained_nuclei.sum()),
        '12w_retained_by_mouse_population': qc_counts.to_dict(orient='records'),
        'gene_selection_population': {
            'genotype': 'WT only', 'ages_weeks': list(DISCOVERY_AGES),
            'n_WT_discovery_nuclei_used': int(WT_DISCOVERY.sum()),
            'n_KI_discovery_nuclei_used': 0, 'n_age12_nuclei_used': 0,
            'excluded_KI_discovery_nuclei_available': int((selected.genotype.eq('KI') & selected.age_weeks.isin(DISCOVERY_AGES)).sum()),
            'excluded_age12_output_nuclei_available': int(selected.age_weeks.eq(12).sum()),
        },
        'T5_ratio_reference': {
            'script_path': T5_RATIO_SCRIPT_REL,
            'script_sha256': T5_RATIO_SCRIPT_SHA256,
            'ratio_function': 'main().ratio_of, s12_t5_decomp_r1.py lines 214-215',
            'pooled_mean_rule': 'WT pooled means over mouse x population units with >=20 WT discovery nuclei, lines 184-192',
        },
    }


def run_real(args: argparse.Namespace) -> dict:
    # G construction has its own gate scope; the imported shared helper bytes are
    # also pinned because this entry point depends on its gate/provenance code.
    expected_scripts = {
        Path(__file__).name: sha(Path(__file__)),
        Path(shared.__file__).name: shared.sha(Path(shared.__file__)),
    }
    gate = shared.check_go(args.go_file, expected_scope='G_builder', expected_scripts=expected_scripts)
    if args.out.exists():
        raise FileExistsError(f'refusing to overwrite existing output: {args.out}')
    env_record, explicit = shared.validate_env_provenance(capture=True)
    pins = check_input_pins(args.h5ad, args.manifest, args.animal_manifest,
                             args.scoping_record, args.full_record, args.ambient_release,
                             args.t5_ratio_script)
    scoping = json.loads(args.scoping_record.read_text(encoding='utf-8'))
    assert scoping['nucleus_manifest']['sha256'] == NUCLEUS_MANIFEST_SHA256
    assert scoping['nucleus_manifest']['rows_of_the_six_mice'] == 38017
    assert scoping['nucleus_manifest']['rows_of_the_six_with_included_true'] == 0
    assert scoping['h5ad']['sha256'] == H5AD_SHA256
    full_record = json.loads(args.full_record.read_text(encoding='utf-8'))
    assert full_record['pins'][_release_paths.expand('${DATA_ROOT}/data/adata_leiden_0.175_annotated_cell_types_v3.h5ad')] == H5AD_SHA256
    manifest = pd.read_csv(args.manifest, sep='\t', dtype={'animal_id': str}, low_memory=False)
    animals = pd.read_csv(args.animal_manifest, sep='\t', dtype={'animal_id': str}, low_memory=False)
    selected, qc_counts = build_manifest_selection(manifest, animals, scoping)
    ambient_ids = ambient_flagged_ids(args.ambient_release)
    with h5py.File(args.h5ad, 'r') as f:
        var = read_elem(f['var'])
        assert 'gene_ids' in var.columns and var.gene_ids.astype(str).is_unique
        counts_shape = f['layers/counts'].attrs.get('shape')
        if counts_shape is not None:
            assert int(counts_shape[1]) == len(var)
    started = time.time()
    result = _build_outputs(args.h5ad, selected, var, ambient_ids, args.out, qc_counts,
                            expected_manifest_obs=set(manifest.obs_name.astype(str)))
    (args.out / 'CONDA_EXPLICIT.txt').write_bytes(explicit)
    env_record['no_environment_changes'] = True
    record = dict(mode='real', scope='all-gene G(a)/G(b) builder', gate_scope='G_builder',
                   approval_gate=gate,
                   script_sha256=expected_scripts, go_file_sha256=sha(args.go_file),
                  input_sha256=pins, nucleus_qc_policy=TWELVE_QC_POLICY,
                   ambient_flag_policy='remove S08/S11 potential_ambient_suspect OR ambient_flagged_hit for the matching population; release table hash-pinned',
                   G_definitions={
                       'G_a': 'genes detected (raw count > 0) in >=1% of WT discovery nuclei of that population',
                       'G_b': 'all-gene ratio >=1 using pooled WT-discovery linear CP10k means; own mean >=1 CP10k and detection >=20% in own WT discovery nuclei; remove T5/S08 ambient-flagged genes',
                       'own_ratio': 'T5 main().ratio_of semantics; pooled WT-discovery means in mouse x population units with >=20 WT discovery nuclei; maximum over other annotated populations',
                       'selection_scope': 'WT discovery nuclei at ages 5,18,24,30 only; no KI or 12 W nucleus enters G(a) or G(b) feature selection',
                       'T5_reference': {'script_path': T5_RATIO_SCRIPT_REL, 'sha256': T5_RATIO_SCRIPT_SHA256,
                                        'ratio_function': 'main().ratio_of lines 214-215; pooled mean convention lines 184-192'},
                       'normalization': 'log1p(10000 * raw gene counts / full-transcriptome raw UMI total)',
                       'fit_population': 'WT discovery ages 5, 18, 24, 30 only; age 12 is output-only'},
                  environment=env_record, result=result, elapsed_seconds=time.time() - started)
    (args.out / 'G_RECORD.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    files = [p for p in sorted(args.out.rglob('*')) if p.is_file() and p.name != 'SHA256SUMS.txt']
    (args.out / 'SHA256SUMS.txt').write_text(''.join(
        f'{sha(p)}  {p.relative_to(args.out).as_posix()}\n' for p in files), encoding='utf-8')
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--real', action='store_true', required=True)
    parser.add_argument('--go-file', type=Path, required=True)
    parser.add_argument('--h5ad', type=Path, default=Path(_release_paths.expand('${DATA_ROOT}/adata_annotated_v3.h5ad')))
    parser.add_argument('--manifest', type=Path, default=Path(_release_paths.expand('${OUTPUT_ROOT}/input_sca1_audit_run01/nucleus_manifest.tsv.gz')))
    parser.add_argument('--animal-manifest', type=Path, default=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'))
    parser.add_argument('--scoping-record', type=Path, default=ROOT / SCOPING_REL)
    parser.add_argument('--full-record', type=Path, default=ROOT / TWELVE_FULL_RECORD_REL)
    parser.add_argument('--ambient-release', type=Path, default=ROOT / AMBIENT_RELEASE_REL)
    parser.add_argument('--t5-ratio-script', type=Path, default=ROOT / T5_RATIO_SCRIPT_REL)
    parser.add_argument('--out', type=Path, default=Path(_release_paths.expand('${WORK_ROOT}/phase4_state_space_G')))
    args = parser.parse_args()
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        raise SystemExit('refusing: launch this entry point with python -I -B')
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('refusing: validation relies on Python assertions')
    envlib = str(shared.ENV / 'lib')
    if envlib not in os.environ.get('LD_LIBRARY_PATH', '').split(':'):
        raise SystemExit('refusing: launch with LD_LIBRARY_PATH=' + envlib + ':$LD_LIBRARY_PATH')
    try:
        result = run_real(args)
        print(json.dumps({'completed': True, 'result': result['result'], 'output': str(args.out)}))
    except BaseException:
        if args.out.exists():
            (args.out / 'FAILED.json').write_text(json.dumps({'traceback': traceback.format_exc()}, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
