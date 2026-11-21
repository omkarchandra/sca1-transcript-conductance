"""Gene-set P-space builder for the Phase 4 shared matrix interface.

Real-data execution is single-shot and byte-gated.  Synthetic tests live in
phase4_state_space_P_fixtures_v05g.py and never open the source H5AD.
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
import re
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Optional, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
from scipy.stats import spearmanr

try:
    from anndata.io import read_elem
except ModuleNotFoundError:
    from anndata._io.specs import read_elem

HERE = Path(__file__).resolve().parent
ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
WORKSPACE = Path(_release_paths.expand('${WORK_ROOT}/state_space'))
POPULATIONS = ('Purkinje cell', 'Bergmann glia')
POP_SLUG = {p: p.replace(' ', '_') for p in POPULATIONS}
P_INTERFACE_TRIPLET_FILES = tuple(
    filename
    for population in POPULATIONS
    for filename in (
        f"state_P_{POP_SLUG[population]}.npy",
        f"state_P_{POP_SLUG[population]}_nuclei.tsv",
        f"state_P_{POP_SLUG[population]}_features.tsv",
    )
)
GRANULE = 'Granule cell'
T6B_ALLPOPS = ('Astrocyte', 'Bergmann glia', 'Endothelial cell', 'GABAergic interneuron 1',
               'GABAergic interneuron 2', 'GABAergic interneuron 3', 'Granule cell', 'Microglia',
               'Oligodendrocyte', 'Oligodendrocyte progenitor cell', 'Pericyte', 'Purkinje cell',
               'Unipolar brush cell')
DISCOVERY_AGES = (5, 18, 24, 30)
READOUT_AGES = (5, 12, 18, 24, 30)
TWELVE_MICE = ('2061', '2062', '22018', '2063', '2065', '22019')
TWELVE_POPULATION_COUNTS = {'Purkinje cell': 203, 'Bergmann glia': 941}
EXPECTED_MICE_PER_POP = 12
H5AD_SHA256 = '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'
NUCLEUS_MANIFEST_SHA256 = 'c3e7d228e97c4991b7d288f46e210308b3ace0cb839aac386f3d578730d08352'
ANIMAL_MANIFEST_SHA256 = '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'
SCOPING_SHA256 = 'a7131e744951ec404f15ba003542a80ad793df814ea17c58f6290a05b8ec7447'
GMT_SHA256 = '712fd1c0f9dcec85b45da2aeb75a355d8f3c22542f0b1a77cd4520d89d2b222d'
BASE_ENV_EXPLICIT_SHA256 = '83a1d457bc19a6e80ad7613bc16ba83a01b8f92a73ab992e8185161950d703a0'
UCELL_VERSION = '2.6.2'
DECOUPLER_VERSION = '1.9.0'
MAX_RANKS = (1000, 1500, 3000)
PRIMARY_MAX_RANK = 1500
MIN_MATCHED_GENES = 15
MAX_MATCHED_GENES = 200
N_RANDOM_SETS = 200
EXPRESSION_BINS = 20
RANDOM_BATCH_FEATURES = 128
CHUNK_ROWS = 256
DEFAULT_UCELL_PARALLEL_WORKERS = 4
DEFAULT_RANDOM_CHUNK_ROWS = 256
SEED = 20261008
ENRICHMENT_SHARE_CUTOFF = 0.95
MIN_ENRICHED_MICE = 11
CORRELATION_FLAG = 0.4
QC4_DOWNSAMPLE_R_FLAG = 0.8
QC5_MAXRANK_JACCARD_FLAG = 0.7
ENV = WORKSPACE / 'env/single_cell_ucell'
CONDA = Path(_release_paths.expand('${CONDA_BIN}'))
H5AD = Path(_release_paths.expand('${DATA_ROOT}/adata_annotated_v3.h5ad'))
MANIFEST = Path(_release_paths.expand('${OUTPUT_ROOT}/input_sca1_audit_run01/nucleus_manifest.tsv.gz'))
ANIMAL_MANIFEST = Path(_release_paths.expand('${OUTPUT_ROOT}/input_sca1_audit_run01/animal_manifest.tsv'))
GMT = WORKSPACE / _release_paths.expand('${DATA_ROOT}/data/msigdb_native_mouse_2025.1.Mm.gmt')
DEFAULT_OUT = WORKSPACE / 'interface/P_v05g_01'
SCOPING = Path(_release_paths.expand('${WORK_ROOT}/SCOPING.json'))
UCELL_HELPER = HERE / 'phase4_state_space_P_ucell_v05g.R'
T6B_SCRIPT = HERE / 's12_t6b_leakin_r2.py'
T6B_SCRIPT_SHA256 = 'b0d57f46604e8d5cef8689b4080a364241d0a1c55facea7688b93e22e453ae72'
T6B_PRIMARY_RATIO = 4.0
T6B_GRANULE_QC6_RATIO = 3.0
T6B_MIN_MEAN_CP10K = 1.0
T6B_MIN_PREVALENCE = 0.20
T6B_TOP_N = 50
T6B_ACTIVE_UNIT_FLOOR = 20
QC6_RHO_CUTOFF = 0.4
QC6_LOO_DRIVER_TOP_FRACTION = 0.10
QC6_LOO_BASE_SET_BATCH = 64

_OFFICIAL_REMOTE = WORKSPACE / _release_paths.expand('${DATA_ROOT}/data/official_msigdb2025.1.Mm')
_OFFICIAL_LOCAL = ROOT / _release_paths.expand('${DATA_ROOT}/analysis_v78/msigdb2025.1.Mm')
_OFFICIAL_ROOT = _OFFICIAL_REMOTE if _OFFICIAL_REMOTE.is_dir() else _OFFICIAL_LOCAL
OFFICIAL_GMTS = {
    'MH': _OFFICIAL_ROOT / 'mh.all.v2025.1.Mm.symbols.gmt',
    'M2_CP_REACTOME': _OFFICIAL_ROOT / 'm2.cp.reactome.v2025.1.Mm.symbols.gmt',
    'M2_CP_WIKIPATHWAYS': _OFFICIAL_ROOT / 'm2.cp.wikipathways.v2025.1.Mm.symbols.gmt',
    'M5_GO_BP': _OFFICIAL_ROOT / 'm5.go.bp.v2025.1.Mm.symbols.gmt',
}


def sha256(path: Union[Path, str]) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def tree_file_size_bytes(directory: Path) -> int:
    return sum(path.stat().st_size for path in directory.rglob('*') if path.is_file())


def write_p_interface_checksum_manifest(interface_dir: Path) -> Path:
    """Bind exactly the six P matrix/nuclei/features files consumed by shared A."""
    missing = [name for name in P_INTERFACE_TRIPLET_FILES if not (interface_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(f'P interface triplet is incomplete; missing: {missing}')
    manifest = interface_dir / 'SHA256SUMS.txt'
    manifest.write_text(''.join(
        f'{sha256(interface_dir / name)}  {name}\n'
        for name in P_INTERFACE_TRIPLET_FILES
    ), encoding='utf-8')
    return manifest


def script_hashes() -> dict[str, str]:
    return {p.name: sha256(p) for p in (Path(__file__), UCELL_HELPER)}


def check_go(path: Path) -> dict[str, Any]:
    gate = json.loads(path.read_text(encoding='utf-8'))
    expected = script_hashes()
    assert gate.get('byte_clearance') == 'yes', 'byte review is absent'
    assert gate.get('author_written_yes') == 'yes', '<AUTHOR> literal written yes is absent'
    assert gate.get('script_sha256') == expected, ('go-file does not bind these exact P bytes', gate.get('script_sha256'), expected)
    assert gate.get('relayed_by') and gate.get('relayed_at_cdt'), 'approval relay provenance is absent'
    assert gate.get('scope') == 'gene_set_space_P', 'go-file scope is not gene_set_space_P'
    assert 'g_interface_sha256sums' not in gate, 'P gate must not bind or depend on the parallel G interface'
    return gate


def parse_gmt(path: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    with path.open(encoding='utf-8') as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 3:
                raise ValueError(f'{path}:{line_number}: GMT row has fewer than 3 columns')
            name = fields[0]
            if name in result:
                raise ValueError(f'{path}:{line_number}: duplicate set name {name}')
            result[name] = sorted(set(g for g in fields[2:] if g))
    return result


def compare_official_gmts(combined: Path, official: dict[str, Path]) -> dict[str, Any]:
    """Compare exact native-mouse set names and member sets; never repair differences."""
    native = parse_gmt(combined)
    details: dict[str, Any] = {}
    for collection, path in official.items():
        expected = parse_gmt(path)
        prefix = collection + '__'
        found = {name[len(prefix):]: genes for name, genes in native.items() if name.startswith(prefix)}
        missing_names = sorted(set(expected) - set(found))
        extra_names = sorted(set(found) - set(expected))
        member_differences = []
        for name in sorted(set(expected) & set(found)):
            missing = sorted(set(expected[name]) - set(found[name]))
            extra = sorted(set(found[name]) - set(expected[name]))
            if missing or extra:
                member_differences.append(dict(set_name=name, missing_members=missing, extra_members=extra))
        details[collection] = dict(native_set_count=len(found), official_set_count=len(expected),
                                   official_gmt_sha256=sha256(path),
                                   missing_set_names=missing_names, extra_set_names=extra_names,
                                   member_difference_count=len(member_differences),
                                   member_differences=member_differences)
    all_equal = all(not row['missing_set_names'] and not row['extra_set_names']
                    and row['member_difference_count'] == 0 for row in details.values())
    return dict(combined_gmt_sha256=sha256(combined), exact_names_and_members_equal=all_equal,
                collections=details)


def build_matched_sets(gmt: dict[str, list[str]], gene_names: np.ndarray) -> tuple[list[dict], pd.DataFrame]:
    gene_set = set(map(str, gene_names))
    records = []
    for full_name, members in sorted(gmt.items()):
        prefix, sep, name = full_name.partition('__')
        if not sep:
            raise ValueError(f'GMT name lacks required collection prefix: {full_name}')
        matched = sorted(gene_set.intersection(members))
        if MIN_MATCHED_GENES <= len(matched) <= MAX_MATCHED_GENES:
            records.append(dict(feature_id=f'P{len(records):06d}', set_name=full_name,
                                collection=prefix, unprefixed_set_name=name,
                                n_native_genes=len(members), n_matched_genes=len(matched),
                                matched_genes=matched))
    feature_table = pd.DataFrame([
        dict(feature_index=i, feature_id=row['feature_id'], stable_id=row['set_name'], set_name=row['set_name'],
             collection=row['collection'], n_native_genes=row['n_native_genes'],
             n_matched_genes=row['n_matched_genes'], matched_gene_symbols=';'.join(row['matched_genes']))
        for i, row in enumerate(records)
    ])
    if feature_table.empty:
        raise ValueError('no MSigDB sets meet the frozen 15–200 matched-gene rule')
    if feature_table.stable_id.duplicated().any():
        raise ValueError('matched feature IDs are not unique')
    return records, feature_table


def expression_bin_ids(mean_cp10k: np.ndarray, gene_names: np.ndarray, n_bins: int = EXPRESSION_BINS) -> np.ndarray:
    mean_cp10k = np.asarray(mean_cp10k, dtype=np.float64)
    assert len(mean_cp10k) == len(gene_names) and np.isfinite(mean_cp10k).all()
    # Equal-count bins with ties broken by gene symbol, so zero-heavy genes are stable.
    order = np.lexsort((np.asarray(gene_names, dtype=str), mean_cp10k))
    bins = np.empty(len(order), dtype=np.int16)
    bins[order] = np.minimum((np.arange(len(order)) * n_bins) // max(1, len(order)), n_bins - 1)
    return bins


def _stable_seed(seed: int, population: str, feature_name: str) -> int:
    digest = hashlib.sha256(f'{seed}|{population}|{feature_name}'.encode()).digest()
    return int.from_bytes(digest[:8], 'little', signed=False)


def make_random_sets(records: list[dict], gene_names: np.ndarray, bins: np.ndarray,
                     population: str, seed: int = SEED, n_random: int = N_RANDOM_SETS) -> list[list[tuple[str, list[str]]]]:
    """Make 200 unique, set-size- and expression-bin-matched controls per set."""
    bin_members = {b: np.flatnonzero(bins == b) for b in np.unique(bins)}
    gene_to_i = {str(g): i for i, g in enumerate(gene_names)}
    output: list[list[tuple[str, list[str]]]] = []
    for row in records:
        member_idx = np.array([gene_to_i[g] for g in row['matched_genes']], dtype=np.int64)
        counts = np.bincount(bins[member_idx], minlength=int(bins.max()) + 1)
        excluded = set(map(int, member_idx))
        available_by_bin = {}
        for b, needed in enumerate(counts):
            if needed:
                pool = bin_members.get(b, np.array([], dtype=np.int64))
                available = pool[~np.isin(pool, np.fromiter(excluded, dtype=np.int64))]
                if len(available) < int(needed):
                    raise ValueError(f'{population}/{row["set_name"]}: expression bin {b} has too few non-set genes')
                available_by_bin[b] = available
        rng = np.random.default_rng(_stable_seed(seed, population, row['set_name']))
        batch = []
        for rep in range(n_random):
            chosen: list[int] = []
            for b, needed in enumerate(counts):
                if not needed:
                    continue
                chosen.extend(map(int, rng.choice(available_by_bin[b], size=int(needed), replace=False)))
            if len(chosen) != len(member_idx) or len(set(chosen)) != len(chosen):
                raise AssertionError('random-set size/uniqueness mismatch')
            batch.append((f'{row["feature_id"]}_R{rep:03d}', sorted(str(gene_names[i]) for i in chosen)))
        output.append(batch)
    return output


def write_random_batches(random_sets: list[list[tuple[str, list[str]]]], out_dir: Path,
                         batch_features: int = RANDOM_BATCH_FEATURES) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=False)
    paths = []
    for first in range(0, len(random_sets), batch_features):
        path = out_dir / f'batch_{first // batch_features:05d}.gmt'
        with path.open('w', encoding='utf-8', newline='\n') as stream:
            for group in random_sets[first:first + batch_features]:
                for random_id, genes in group:
                    stream.write('\t'.join([random_id, 'matched-random-control', *genes]) + '\n')
        paths.append(path)
    return paths


def write_random_batches_for_population(records: list[dict], gene_names: np.ndarray, bins: np.ndarray,
                                       population: str, out_dir: Path,
                                       batch_features: int = RANDOM_BATCH_FEATURES) -> list[Path]:
    """Generate controls a batch at a time so the full 200×set design stays off RAM."""
    out_dir.mkdir(parents=True, exist_ok=False)
    paths = []
    for first in range(0, len(records), batch_features):
        batch = make_random_sets(records[first:first + batch_features], gene_names, bins, population)
        path = out_dir / f'batch_{first // batch_features:05d}.gmt'
        with path.open('w', encoding='utf-8', newline='\n') as stream:
            for group in batch:
                for random_id, genes in group:
                    stream.write('\t'.join([random_id, 'matched-random-control', *genes]) + '\n')
        paths.append(path)
    return paths


def summarize_enrichment(mouse_median_shares_by_rank: dict[int, np.ndarray], ulm_by_mouse: np.ndarray,
                         mouse_ids: list[str], feature_ids: list[str]) -> tuple[pd.DataFrame, dict[int, set[str]]]:
    ulm_by_mouse = np.asarray(ulm_by_mouse, dtype=np.float64)
    if not np.isfinite(ulm_by_mouse).all():
        raise ValueError('ULM values must be finite before enrichment')
    rows = []
    enriched: dict[int, set[str]] = {}
    for rank, mouse_median_shares in mouse_median_shares_by_rank.items():
        assert mouse_median_shares.shape == ulm_by_mouse.shape == (len(mouse_ids), len(feature_ids))
        rank_hits = set()
        pass_by_feature = {}
        for j, feature_id in enumerate(feature_ids):
            share_mice = mouse_median_shares[:, j] >= ENRICHMENT_SHARE_CUTOFF
            ulm_mice = ulm_by_mouse[:, j] > 0
            same_mice = share_mice & ulm_mice
            pass_feature = int(same_mice.sum()) >= MIN_ENRICHED_MICE
            pass_by_feature[feature_id] = pass_feature
            if pass_feature:
                rank_hits.add(feature_id)
            for i, mouse in enumerate(mouse_ids):
                rows.append(dict(maxRank=rank, feature_id=feature_id, mouse_id=mouse,
                                 median_nucleus_share_of_random_sets_beaten=float(mouse_median_shares[i, j]),
                                 share_pass=bool(share_mice[i]), ULM_estimate=float(ulm_by_mouse[i, j]),
                                 ULM_positive=bool(ulm_mice[i]), same_mouse_pass=bool(same_mice[i]),
                                 enriched_set_pass=bool(pass_by_feature[feature_id])))
        enriched[rank] = rank_hits
    return pd.DataFrame(rows), enriched


def enrichment_from_nucleus_shares(nucleus_shares: np.ndarray, meta: pd.DataFrame,
                                   mouse_ids: list[str]) -> np.ndarray:
    """Median per-nucleus random-control share for each mouse × feature."""
    nucleus_shares = np.asarray(nucleus_shares, dtype=np.float64)
    assert nucleus_shares.ndim == 2 and nucleus_shares.shape[0] == len(meta)
    out = np.empty((len(mouse_ids), nucleus_shares.shape[1]), dtype=np.float64)
    labels = meta.mouse_id.astype(str).to_numpy()
    wt = meta.genotype.astype(str).str.upper().eq('WT').to_numpy()
    discovery = meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()
    for i, mouse in enumerate(mouse_ids):
        rows = (labels == str(mouse)) & wt & discovery
        if not rows.any():
            raise ValueError(f'no WT discovery nuclei for mouse {mouse}')
        out[i] = np.median(nucleus_shares[rows], axis=0)
    return out


def qc_hard_gates(scores: np.ndarray, iqr_zero_fraction: Optional[float] = None,
                  sets_ge_099_per_nucleus: Optional[np.ndarray] = None) -> dict[str, Any]:
    """QC1/QC2 are the only hard gates; QC3–QC6 remain reporting flags."""
    matrix = np.asarray(scores, dtype=np.float64)
    if matrix.ndim != 2 or not matrix.size:
        raise ValueError('QC gate input must be a finite nonempty nuclei × feature matrix')
    if sets_ge_099_per_nucleus is None:
        counts = []
        for first in range(0, len(matrix), CHUNK_ROWS):
            block = np.asarray(matrix[first:first + CHUNK_ROWS], dtype=np.float64)
            if not np.isfinite(block).all():
                raise ValueError('QC gate input has non-finite scores')
            counts.append(np.sum(block >= 0.99, axis=1))
        n_near_ceiling_per_nucleus = np.concatenate(counts)
    else:
        n_near_ceiling_per_nucleus = np.asarray(sets_ge_099_per_nucleus, dtype=np.int64)
        assert n_near_ceiling_per_nucleus.shape == (len(matrix),)
    if iqr_zero_fraction is None:
        quartiles = np.quantile(matrix, [0.25, 0.75], axis=1)
        iqr_zero_fraction = float(np.mean((quartiles[1] - quartiles[0]) == 0.0))
    qc1_median = float(np.median(n_near_ceiling_per_nucleus))
    failures = []
    if qc1_median > 0:
        failures.append('QC1: median number of sets per nucleus scoring >=0.99 is greater than zero')
    if float(iqr_zero_fraction) > 0.01:
        failures.append('QC2: more than 1% of nuclei have within-nucleus score IQR equal to zero')
    return dict(QC1_median_sets_ge_0_99=qc1_median,
                QC2_fraction_nuclei_IQR_zero=float(iqr_zero_fraction),
                QC1_status='STOP' if qc1_median > 0 else 'PASS',
                QC2_status='STOP' if float(iqr_zero_fraction) > 0.01 else 'PASS',
                hard_gate_status='STOP' if failures else 'PASS', hard_gate_reasons=failures)


def group_redundant_sets(enriched_ids: set[str], records: list[dict], wt_discovery_scores: np.ndarray,
                         feature_ids: list[str], jaccard_cutoff: float = 0.5,
                         correlation_cutoff: float = 0.9) -> tuple[pd.DataFrame, pd.DataFrame]:
    ordered = sorted(enriched_ids)
    if not ordered:
        return pd.DataFrame(columns=['group_id', 'representative_feature_id', 'member_feature_id', 'n_matched_genes', 'jaccard_to_representative', 'WT_discovery_nuclei_score_correlation_to_representative']), pd.DataFrame(columns=['group_id', 'representative_feature_id', 'n_members', 'n_matched_genes'])
    rec_by_id = {row['feature_id']: row for row in records}
    col_by_id = {name: i for i, name in enumerate(feature_ids)}
    parent = list(range(len(ordered)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    sets = [set(rec_by_id[x]['matched_genes']) for x in ordered]
    vectors = np.asarray(wt_discovery_scores[:, [col_by_id[x] for x in ordered]], dtype=np.float64)
    centered = vectors - vectors.mean(axis=0, keepdims=True)
    norms = np.sqrt(np.sum(centered * centered, axis=0))
    for i in range(len(ordered)):
        for j in range(i + 1, len(ordered)):
            union_pair = False
            inter = len(sets[i] & sets[j])
            jac = inter / len(sets[i] | sets[j]) if sets[i] | sets[j] else 0.0
            if jac >= jaccard_cutoff:
                union_pair = True
            elif norms[i] > 0 and norms[j] > 0:
                rho = float(np.dot(centered[:, i], centered[:, j]) / (norms[i] * norms[j]))
                union_pair = rho >= correlation_cutoff
            if union_pair:
                union(i, j)
    groups: dict[int, list[int]] = {}
    for i in range(len(ordered)):
        groups.setdefault(find(i), []).append(i)
    group_rows, representative_rows = [], []
    sorted_groups = sorted(groups.values(), key=lambda ix: min(ordered[i] for i in ix))
    for group_number, members in enumerate(sorted_groups, start=1):
        representative = min((ordered[i] for i in members), key=lambda name: (len(rec_by_id[name]['matched_genes']), name))
        rep_i = ordered.index(representative)
        group_id = f'G{group_number:05d}'
        representative_rows.append(dict(group_id=group_id, representative_feature_id=representative,
                                         n_members=len(members), n_matched_genes=len(rec_by_id[representative]['matched_genes'])))
        for i in sorted(members, key=lambda ix: ordered[ix]):
            name = ordered[i]
            inter = len(sets[rep_i] & sets[i])
            jac = inter / len(sets[rep_i] | sets[i]) if sets[rep_i] | sets[i] else 1.0
            if norms[rep_i] > 0 and norms[i] > 0:
                corr = float(np.dot(centered[:, rep_i], centered[:, i]) / (norms[rep_i] * norms[i]))
            else:
                corr = np.nan
            group_rows.append(dict(group_id=group_id, representative_feature_id=representative,
                                   member_feature_id=name, n_matched_genes=len(sets[i]),
                                   jaccard_to_representative=jac,
                                    WT_discovery_nuclei_score_correlation_to_representative=corr))
    return pd.DataFrame(group_rows), pd.DataFrame(representative_rows)


def downsample_mouse_means(scores: np.ndarray, meta: pd.DataFrame, mouse_ids: list[str],
                           seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """Per-mouse 50% WT-discovery resample and full-vs-resampled set means."""
    rng = np.random.default_rng(seed)
    wt_mask = meta.genotype.astype(str).str.upper().eq('WT').to_numpy() & meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()
    wt_meta = meta.loc[wt_mask].reset_index(drop=True)
    wt_scores = np.asarray(scores[wt_mask], dtype=np.float64)
    full_means = np.empty((len(mouse_ids), wt_scores.shape[1]), dtype=np.float64)
    down_means = np.empty_like(full_means)
    for mouse_i, mouse in enumerate(mouse_ids):
        local = np.flatnonzero(wt_meta.mouse_id.astype(str).to_numpy() == mouse)
        if not len(local):
            raise ValueError(f'WT discovery mouse has no nuclei for downsample check: {mouse}')
        n_take = max(1, int(np.ceil(len(local) * 0.5)))
        chosen = np.sort(rng.choice(local, size=n_take, replace=False))
        full_means[mouse_i] = wt_scores[local].mean(axis=0)
        down_means[mouse_i] = wt_scores[chosen].mean(axis=0)
    return full_means, down_means


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')


def _strict_manifest_boolean(series: pd.Series, column: str) -> np.ndarray:
    labels = series.astype(str).str.strip().str.upper()
    unexpected = sorted(set(labels) - {'TRUE', 'FALSE'})
    if unexpected:
        raise ValueError(f'{column} must contain only TRUE/FALSE labels, found {unexpected}')
    return labels.eq('TRUE').to_numpy()


def build_p_nucleus_metadata(h5ad: Path, manifest: pd.DataFrame, animals: pd.DataFrame,
                             scoping: dict[str, Any],
                             strict_expected_counts: bool = True) -> dict[str, pd.DataFrame]:
    """Independently reproduce the nucleus selection and sort each output by H5AD row."""
    required_manifest = {'obs_name', 'animal_id', 'cell_type_raw', 'included', 'exclusion_reason',
                         'age_weeks_obs', 'genotype_final'}
    required_animals = {'animal_id', 'age_weeks', 'genotype', 'sex'}
    if not required_manifest <= set(manifest.columns):
        raise ValueError(f'nucleus manifest lacks columns: {sorted(required_manifest - set(manifest.columns))}')
    if not required_animals <= set(animals.columns):
        raise ValueError(f'animal manifest lacks columns: {sorted(required_animals - set(animals.columns))}')
    if manifest.obs_name.astype(str).duplicated().any():
        raise ValueError('nucleus manifest contains duplicate obs_name values')
    animals = animals.copy()
    animals['animal_id'] = animals.animal_id.astype(str)
    if animals.animal_id.duplicated().any():
        raise ValueError('animal manifest contains duplicate animal_id values')
    manifest = manifest.copy()
    manifest['animal_id'] = manifest.animal_id.astype(str)
    manifest['age_weeks_obs'] = pd.to_numeric(manifest.age_weeks_obs, errors='raise').astype(int)
    included = _strict_manifest_boolean(manifest.included, 'included')
    expected_six = scoping.get('six_mice', {})
    if set(expected_six) != set(TWELVE_MICE):
        raise ValueError('pinned scoping record does not describe the six locked 12-week mice')
    if strict_expected_counts:
        if scoping.get('nucleus_manifest', {}).get('sha256') != NUCLEUS_MANIFEST_SHA256:
            raise ValueError('SCOPING.json is not bound to the pinned nucleus manifest')
        if scoping.get('h5ad', {}).get('sha256') != H5AD_SHA256:
            raise ValueError('SCOPING.json is not bound to the pinned source H5AD')
        if scoping.get('nucleus_manifest', {}).get('rows_of_the_six_mice') != 38017:
            raise ValueError('SCOPING.json six-mouse manifest total changed')
        if scoping.get('nucleus_manifest', {}).get('rows_of_the_six_with_included_true') != 0:
            raise ValueError('SCOPING.json no longer records all six mice as age-excluded')
    if not set(TWELVE_MICE) <= set(animals.animal_id.astype(str)):
        raise ValueError('animal manifest is missing one or more scoped 12-week mice')

    twelve_all = manifest.animal_id.isin(TWELVE_MICE).to_numpy()
    twelve_rows = manifest.loc[twelve_all]
    if strict_expected_counts and len(twelve_rows) != 38017:
        raise ValueError(f'expected 38,017 manifest rows for the six 12-week mice; found {len(twelve_rows)}')
    if included[twelve_all].any():
        raise ValueError('a scoped 12-week row is already marked included; expected age-only exclusion')
    if not twelve_rows.exclusion_reason.astype(str).eq('age_weeks in design.excluded_ages').all():
        raise ValueError('a scoped 12-week row has an exclusion reason other than the design age exclusion')
    if not twelve_rows.age_weeks_obs.eq(12).all():
        raise ValueError('a scoped 12-week mouse has a manifest age other than 12 weeks')
    if not twelve_rows.genotype_final.astype(str).str.upper().isin({'WT', 'KI'}).all():
        raise ValueError('a scoped 12-week row has an unexpected genotype')

    animal_lookup = animals.set_index('animal_id', verify_integrity=True)
    for animal_id, expected in expected_six.items():
        rows = twelve_rows.loc[twelve_rows.animal_id.eq(str(animal_id))]
        animal = animal_lookup.loc[str(animal_id)]
        if len(rows) != int(expected['n_nuclei_table_s1']):
            raise ValueError(f'{animal_id}: manifest row count differs from the scoping record')
        if int(animal.age_weeks) != 12 or str(animal.genotype).upper() != str(expected['genotype']).upper():
            raise ValueError(f'{animal_id}: animal-manifest age/genotype differs from the scoping record')
        if str(animal.sex).upper() != str(expected['sex']).upper():
            raise ValueError(f'{animal_id}: animal-manifest sex differs from the scoping record')

    discovery = included & manifest.age_weeks_obs.isin(DISCOVERY_AGES).to_numpy()
    target = manifest.cell_type_raw.astype(str).isin(POPULATIONS).to_numpy()
    twelve_target = twelve_all & target
    selection = (discovery & target) | twelve_target
    selected = manifest.loc[selection].copy().reset_index(drop=True)
    selected = selected.merge(
        animals[['animal_id', 'age_weeks', 'genotype', 'sex']], on='animal_id', how='left',
        validate='many_to_one', suffixes=('', '_animal'))
    if selected[['age_weeks', 'genotype', 'sex']].isna().any().any():
        raise ValueError('selected P nuclei could not be joined to the animal manifest')
    if not selected.age_weeks.astype(int).eq(selected.age_weeks_obs.astype(int)).all():
        raise ValueError('selected nucleus age differs between the animal and nucleus manifests')
    if not selected.genotype.astype(str).str.upper().eq(selected.genotype_final.astype(str).str.upper()).all():
        raise ValueError('selected nucleus genotype differs between the animal and nucleus manifests')
    selected['age_weeks'] = selected.age_weeks.astype(int)
    selected['genotype'] = selected.genotype.astype(str).str.upper()
    selected['sex'] = selected.sex.astype(str).str.upper()
    selected['animal_id'] = selected.animal_id.astype(str)
    if not selected.age_weeks.isin(READOUT_AGES).all() or not selected.genotype.isin({'WT', 'KI'}).all():
        raise ValueError('selected P nuclei violate the shared readout age/genotype contract')

    twelve_target_rows = selected.loc[selected.age_weeks.eq(12)]
    for population in POPULATIONS:
        population_12w = twelve_target_rows.loc[twelve_target_rows.cell_type_raw.eq(population)]
        if strict_expected_counts and len(population_12w) != TWELVE_POPULATION_COUNTS[population]:
            raise ValueError(f'{population}: expected {TWELVE_POPULATION_COUNTS[population]} scoped 12-week nuclei; found {len(population_12w)}')
        retained = population_12w.groupby('animal_id', sort=True).size()
        expected_counts = scoping['per_mouse_by_celltype_counts']
        for animal_id in TWELVE_MICE:
            expected_n = int(expected_counts[animal_id][population]['nuclei'])
            if int(retained.get(animal_id, 0)) != expected_n:
                raise ValueError(f'{animal_id}/{population}: retained 12-week count differs from SCOPING.json')

    with h5py.File(h5ad, 'r') as handle:
        obs = read_elem(handle['obs'])
        obs_names = obs.index.astype(str).to_numpy()
        if not obs.index.is_unique:
            raise ValueError('H5AD obs names are not unique')
    source_rows = pd.Index(obs_names).get_indexer(selected.obs_name.astype(str))
    if (source_rows < 0).any() or len(np.unique(source_rows)) != len(source_rows):
        raise ValueError('selected P nucleus IDs do not map uniquely to the H5AD obs axis')
    selected['source_row'] = source_rows
    columns = ['nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex']
    selected = selected.rename(columns={'obs_name': 'nucleus_id', 'animal_id': 'mouse_id'})
    result = {}
    for population in POPULATIONS:
        frame = selected.loc[selected.cell_type_raw.astype(str).eq(population)].copy()
        frame = frame.sort_values('source_row', kind='mergesort')
        metadata = frame[columns].reset_index(drop=True)
        if metadata.nucleus_id.duplicated().any():
            raise ValueError(f'{population}: duplicate P nucleus IDs after selection')
        if set(metadata.age_weeks.astype(int)) != set(READOUT_AGES) or set(metadata.genotype) != {'WT', 'KI'}:
            raise ValueError(f'{population}: P nucleus selection does not cover the shared readout contract')
        result[population] = metadata
    return result


def _read_checksum_manifest(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line_number, line in enumerate(path.read_text(encoding='utf-8').splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r'[0-9a-fA-F]{64}', parts[0]):
            raise ValueError(f'{path}:{line_number}: malformed SHA-256 manifest row')
        rel = parts[1].strip().lstrip('*').removeprefix('./')
        if not rel or rel in checksums:
            raise ValueError(f'{path}:{line_number}: blank or duplicate path in SHA-256 manifest')
        checksums[rel] = parts[0].lower()
    if not checksums:
        raise ValueError(f'{path}: empty SHA-256 manifest')
    return checksums


def _read_t6b_wt_manifest(manifest_path: Path, animal_path: Path) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path, sep='\t', dtype={'animal_id': str}, low_memory=False)
    animals = pd.read_csv(animal_path, sep='\t', dtype={'animal_id': str}, low_memory=False)
    req = {'obs_name', 'animal_id', 'cell_type_raw', 'included'}
    assert req <= set(manifest.columns), req - set(manifest.columns)
    keep = (manifest.cell_type_raw.astype(str).isin(T6B_ALLPOPS)
            & manifest.included.astype(str).str.lower().eq('true'))
    cells = manifest.loc[keep, ['obs_name', 'animal_id', 'cell_type_raw']].copy()
    cells = cells.merge(animals[['animal_id', 'age_weeks', 'genotype', 'sex']],
                        on='animal_id', how='left', validate='many_to_one')
    cells['age_weeks'] = cells.age_weeks.astype(int)
    cells['genotype'] = cells.genotype.astype(str).str.upper()
    cells = cells[cells.age_weeks.isin(DISCOVERY_AGES) & cells.genotype.eq('WT')].copy()
    # T6b defines active mouse × population units as at least 20 retained nuclei.
    unit_sizes = cells.groupby(['animal_id', 'cell_type_raw']).size()
    active = {(str(animal), str(pop)) for (animal, pop), n in unit_sizes.items()
              if int(n) >= T6B_ACTIVE_UNIT_FLOOR}
    cells = cells[[ (str(a), str(p)) in active
                    for a, p in zip(cells.animal_id, cells.cell_type_raw) ]].copy()
    cells = cells.rename(columns={'obs_name': 'nucleus_id', 'animal_id': 'mouse_id'})
    cells['mouse_id'] = cells.mouse_id.astype(str)
    cells['sex'] = cells.sex.astype(str).str.upper()
    cells['log10_total_umi'] = np.nan
    cells['n_detected'] = -1
    assert len(cells) and cells.nucleus_id.is_unique
    missing = set(T6B_ALLPOPS) - set(cells.cell_type_raw.astype(str))
    assert not missing, f'T6b marker reference lacks active WT discovery units for {sorted(missing)}'
    return cells[['nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex', 'cell_type_raw',
                  'log10_total_umi', 'n_detected']].reset_index(drop=True)


def compute_t6b_top_markers(h5ad: Path, wt_cells: pd.DataFrame,
                            chunk_rows: int = 2048,
                            granule_ratio: float = T6B_GRANULE_QC6_RATIO) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Recompute pinned T6b WT markers; only the predeclared QC6 Granule ratio is relaxed."""
    if not T6B_SCRIPT.is_file() or sha256(T6B_SCRIPT) != T6B_SCRIPT_SHA256:
        raise RuntimeError('pinned T6b marker-definition script is absent or its source bytes changed')
    with h5py.File(h5ad, 'r') as handle:
        obs = read_elem(handle['obs'])
        var = read_elem(handle['var'])
        obs_names = obs.index.astype(str).to_numpy()
        symbols = var.index.astype(str).to_numpy()
        stable_ids = var['gene_ids'].astype(str).to_numpy()
        assert pd.Index(stable_ids).is_unique and len(stable_ids) == len(symbols)
        row_index = pd.Index(obs_names).get_indexer(wt_cells.nucleus_id.astype(str))
        assert (row_index >= 0).all() and len(np.unique(row_index)) == len(row_index)
        order = np.argsort(row_index, kind='mergesort')
        wt_cells = wt_cells.iloc[order].reset_index(drop=True)
        row_index = row_index[order]
        source_to_wt = np.full(len(obs_names), -1, dtype=np.int64)
        source_to_wt[row_index] = np.arange(len(row_index), dtype=np.int64)
        pop_index = {p: i for i, p in enumerate(T6B_ALLPOPS)}
        pop_ids = np.array([pop_index[p] for p in wt_cells.cell_type_raw.astype(str)], dtype=np.int16)
        sums = np.zeros((len(T6B_ALLPOPS), len(symbols)), dtype=np.float64)
        detections = np.zeros((len(T6B_ALLPOPS), len(symbols)), dtype=np.int64)
        n_cells = np.bincount(pop_ids, minlength=len(T6B_ALLPOPS)).astype(np.int64)
        counts = handle['layers/counts']
        ptr_ds, data_ds, index_ds = counts['indptr'], counts['data'], counts['indices']
        for start in range(0, len(obs_names), chunk_rows):
            stop = min(len(obs_names), start + chunk_rows)
            out_idx = source_to_wt[start:stop]
            local_rows = np.flatnonzero(out_idx >= 0)
            if not len(local_rows):
                continue
            ptr = ptr_ds[start:stop + 1].astype(np.int64)
            lo, hi = int(ptr[0]), int(ptr[-1])
            block_ptr = (ptr - lo).astype(np.int64)
            values = data_ds[lo:hi].astype(np.float64)
            indices = index_ds[lo:hi].astype(np.int64)
            assert np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.floor(values)).all()
            block = sp.csr_matrix((values, indices, block_ptr), shape=(stop - start, len(symbols)))
            block.sum_duplicates(); block.eliminate_zeros()
            block = block[local_rows].tocsr()
            chosen_out_idx = out_idx[local_rows]
            totals = np.asarray(block.sum(axis=1)).ravel()
            assert (totals > 0).all()
            for p_i in range(len(T6B_ALLPOPS)):
                rows = np.flatnonzero(pop_ids[chosen_out_idx] == p_i)
                if not len(rows):
                    continue
                part = block[rows].tocsr()
                cp = part.multiply((10000.0 / totals[rows])[:, None]).tocsr()
                sums[p_i] += np.asarray(cp.sum(axis=0)).ravel()
                detections[p_i] += np.asarray(part.getnnz(axis=0)).ravel().astype(np.int64)
    means = sums / n_cells[:, None]
    prevalence = detections / n_cells[:, None]
    rows = []
    marker_symbols: dict[str, list[str]] = {}
    for p_i, population in enumerate(T6B_ALLPOPS):
        others = [i for i in range(len(T6B_ALLPOPS)) if i != p_i]
        max_other = np.max(means[others], axis=0)
        own = means[p_i]
        ratio = np.divide(own, max_other, out=np.full_like(own, np.inf), where=max_other > 0)
        threshold = float(granule_ratio) if population == GRANULE else T6B_PRIMARY_RATIO
        eligible = np.flatnonzero((own >= T6B_MIN_MEAN_CP10K)
                                  & (prevalence[p_i] >= T6B_MIN_PREVALENCE)
                                  & (ratio >= threshold))
        ranked = eligible[np.lexsort((eligible, -ratio[eligible]))][:T6B_TOP_N]
        symbols_top = list(dict.fromkeys(str(symbols[i]) for i in ranked))
        marker_symbols[population] = symbols_top
        for rank, feature_i in enumerate(ranked, start=1):
            rows.append(dict(population=population, rank=rank, feature_index=int(feature_i),
                             stable_id=str(stable_ids[feature_i]), symbol=str(symbols[feature_i]),
                             mean_WT_CP10k=float(own[feature_i]),
                             WT_detection_prevalence=float(prevalence[p_i, feature_i]),
                             own_vs_max_other_mean_ratio=float(ratio[feature_i]),
                             threshold=threshold))
    marker_table = pd.DataFrame(rows)
    return marker_table, marker_symbols


def t6b_bergmann_neighbour_symbols(marker_symbols: dict[str, list[str]]) -> list[str]:
    """T6b primary Granule/Purkinje marker union after Bergmann self-marker exclusion."""
    required = {GRANULE, 'Purkinje cell', 'Bergmann glia'}
    if not required <= set(marker_symbols):
        raise ValueError(f'T6b marker table lacks populations: {sorted(required - set(marker_symbols))}')
    neighbour = set(marker_symbols[GRANULE]) | set(marker_symbols['Purkinje cell'])
    neighbour -= set(marker_symbols['Bergmann glia'])
    return sorted(neighbour)


def qc5_jaccard_flag(jaccard: float, cutoff: float = QC5_MAXRANK_JACCARD_FLAG) -> bool:
    return bool(not np.isfinite(jaccard) or float(jaccard) < cutoff)


def _write_mtx_chunk(directory: Path, number: int, counts: sp.csr_matrix,
                     metadata: pd.DataFrame) -> None:
    prefix = directory / f'chunk_{number:06d}'
    mmwrite(str(prefix) + '.mtx', counts.astype(np.int64), field='integer')
    metadata.to_csv(str(prefix) + '.nuclei.tsv', sep='\t', index=False)


def stage_source(h5ad: Path, nuclei_by_pop: dict[str, pd.DataFrame], granule: pd.DataFrame,
                 neighbor_marker_symbols: list[str], work: Path,
                 chunk_rows: int = CHUNK_ROWS,
                 random_chunk_rows: Optional[int] = None) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, dict[str, np.ndarray]], dict[str, dict[str, int]], dict[str, Any]]:
    """Stream pinned CSR counts into observed and WT-random MatrixMarket chunks."""
    if random_chunk_rows is None:
        random_chunk_rows = chunk_rows
    if chunk_rows < 1 or random_chunk_rows < 1:
        raise ValueError('observed and WT-random chunk row counts must be positive')
    work.mkdir(parents=True, exist_ok=False)
    target_work = {p: work / POP_SLUG[p] for p in POPULATIONS}
    for p in POPULATIONS:
        target_work[p].mkdir()
    wt_random_work = {p: work / f'{POP_SLUG[p]}_wt_discovery_random' for p in POPULATIONS}
    for directory in wt_random_work.values():
        directory.mkdir()
    granule_work = work / 'granule_leak_only'
    granule_work.mkdir()

    with h5py.File(h5ad, 'r') as handle:
        obs = read_elem(handle['obs'])
        var = read_elem(handle['var'])
        obs_names = obs.index.astype(str).to_numpy()
        assert obs.index.is_unique and 'gene_ids' in var.columns
        stable = var.gene_ids.astype(str).to_numpy()
        symbols = var.index.astype(str).to_numpy()
        assert len(stable) == len(symbols) and pd.Index(stable).is_unique and np.all(stable != '')
        assert np.all(symbols != '')
        gene_names = np.array(sorted(set(symbols)), dtype=str)
        gene_to_index = {name: i for i, name in enumerate(gene_names)}
        duplicate_symbol_counts = pd.Series(symbols).value_counts()
        duplicate_symbol_names = sorted(duplicate_symbol_counts[duplicate_symbol_counts > 1].index.astype(str))
        n_duplicate_symbol_columns_aggregated = int(len(symbols) - len(gene_names))
        original_to_symbol = np.fromiter((gene_to_index[s] for s in symbols), dtype=np.int64, count=len(symbols))
        aggregate = sp.csr_matrix((np.ones(len(symbols), dtype=np.int8),
                                   (np.arange(len(symbols), dtype=np.int64), original_to_symbol)),
                                  shape=(len(symbols), len(gene_names)))
        counts_group = handle['layers/counts']
        shape = tuple(int(x) for x in counts_group.attrs['shape'])
        assert shape == (len(obs_names), len(symbols))
        ptr_ds, data_ds, index_ds = counts_group['indptr'], counts_group['data'], counts_group['indices']

        metadata_by_label: dict[str, pd.DataFrame] = {p: nuclei_by_pop[p].copy() for p in POPULATIONS}
        metadata_by_label[GRANULE] = granule.copy()
        for frame in metadata_by_label.values():
            frame['log10_total_umi'] = np.full(len(frame), np.nan, dtype=np.float64)
            frame['n_detected'] = np.full(len(frame), -1, dtype=np.int64)
        marker_indices = np.array(sorted({gene_to_index[g] for g in neighbor_marker_symbols if g in gene_to_index}), dtype=np.int64)
        if not len(marker_indices):
            raise ValueError('no T6b-defined Granule/Purkinje marker symbols are present in the P gene axis')
        marker_indicator = np.zeros(len(gene_names), dtype=np.float64)
        marker_indicator[marker_indices] = 1.0
        source_map: dict[str, np.ndarray] = {}
        for label, frame in metadata_by_label.items():
            row_index = pd.Index(obs_names).get_indexer(frame.nucleus_id.astype(str))
            assert (row_index >= 0).all() and len(np.unique(row_index)) == len(row_index)
            if label in POPULATIONS:
                assert np.all(np.diff(row_index) > 0), f'{label}: manifest metadata is not in H5AD source order'
            else:
                order = np.argsort(row_index, kind='mergesort')
                frame = frame.iloc[order].reset_index(drop=True)
                row_index = row_index[order]
            assert np.all(np.diff(row_index) > 0), f'{label}: selected nuclei are not in source order'
            metadata_by_label[label] = frame
            src_to_out = np.full(len(obs_names), -1, dtype=np.int64)
            src_to_out[row_index] = np.arange(len(row_index), dtype=np.int64)
            source_map[label] = src_to_out

        # Raw-count metadata is populated during the source pass before duplicate-symbol aggregation.
        bergmann_neighbour_score = np.full(len(metadata_by_label['Bergmann glia']), np.nan, dtype=np.float64)

        pending_counts = {label: sp.csr_matrix((0, len(gene_names)), dtype=np.float64)
                          for label in (*POPULATIONS, GRANULE)}
        pending_meta = {label: metadata_by_label[label].iloc[:0].copy()
                        for label in (*POPULATIONS, GRANULE)}
        chunk_number = {label: 0 for label in (*POPULATIONS, GRANULE)}
        pending_wt_counts = {p: sp.csr_matrix((0, len(gene_names)), dtype=np.float64)
                             for p in POPULATIONS}
        pending_wt_meta = {p: metadata_by_label[p].iloc[:0].copy() for p in POPULATIONS}
        wt_chunk_number = {p: 0 for p in POPULATIONS}
        wt_random_rows_written = {p: 0 for p in POPULATIONS}
        mean_cp10k_sum = {p: np.zeros(len(gene_names), dtype=np.float64) for p in POPULATIONS}
        wt_cell_count = {p: 0 for p in POPULATIONS}
        mouse_ids = {p: sorted(metadata_by_label[p].loc[
            metadata_by_label[p].genotype.eq('WT') & metadata_by_label[p].age_weeks.isin(DISCOVERY_AGES), 'mouse_id'].astype(str).unique())
            for p in POPULATIONS}
        ulm_sum = {p: {m: np.zeros(len(gene_names), dtype=np.float64) for m in mouse_ids[p]} for p in POPULATIONS}
        ulm_n = {p: {m: 0 for m in mouse_ids[p]} for p in POPULATIONS}

        for start in range(0, len(obs_names), 2048):
            stop = min(len(obs_names), start + 2048)
            ptr = ptr_ds[start:stop + 1].astype(np.int64)
            lo, hi = int(ptr[0]), int(ptr[-1])
            local_ptr = ptr - lo
            values = data_ds[lo:hi].astype(np.float64)
            cols = index_ds[lo:hi].astype(np.int64)
            assert np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.floor(values)).all()
            row_idx = np.repeat(np.arange(stop - start, dtype=np.int64), np.diff(local_ptr))
            raw = sp.csr_matrix((values, (row_idx, cols)), shape=(stop - start, len(symbols)))
            raw.sum_duplicates(); raw.eliminate_zeros()
            raw_unique = (raw @ aggregate).tocsr()
            totals = np.asarray(raw.sum(axis=1)).ravel()
            detected = np.asarray(raw.getnnz(axis=1)).ravel().astype(np.int64)
            assert (totals > 0).all()

            for label in (*POPULATIONS, GRANULE):
                dest = source_map[label][start:stop]
                keep_rows = np.flatnonzero(dest >= 0)
                if not len(keep_rows):
                    continue
                out_idx = dest[keep_rows]
                raw_block = raw_unique[keep_rows].tocsr()
                local_totals = totals[keep_rows]
                local_detected = detected[keep_rows]
                metadata_by_label[label].loc[out_idx, 'log10_total_umi'] = np.log10(local_totals)
                metadata_by_label[label].loc[out_idx, 'n_detected'] = local_detected
                frame = metadata_by_label[label].iloc[out_idx].copy().reset_index(drop=True)
                if label in POPULATIONS:
                    assert np.allclose(np.log10(local_totals), frame.log10_total_umi.to_numpy(dtype=float), rtol=0, atol=1e-12)
                    assert np.array_equal(local_detected, frame.n_detected.to_numpy(dtype=np.int64))
                    discovery = frame.genotype.eq('WT').to_numpy() & frame.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()
                    if discovery.any():
                        cp = raw_block[discovery].multiply((10000.0 / local_totals[discovery])[:, None]).tocsr()
                        mean_cp10k_sum[label] += np.asarray(cp.sum(axis=0)).ravel()
                        wt_cell_count[label] += int(discovery.sum())
                        log_cp = cp.copy(); log_cp.data = np.log1p(log_cp.data)
                        if label == 'Bergmann glia':
                            neighbor_rows = np.flatnonzero(discovery)
                            neighbor_values = np.asarray(log_cp @ marker_indicator).ravel()
                            bergmann_neighbour_score[out_idx[neighbor_rows]] = neighbor_values
                        wt_meta = frame.iloc[np.flatnonzero(discovery)]
                        for mouse in mouse_ids[label]:
                            chosen = np.flatnonzero(wt_meta.mouse_id.astype(str).to_numpy() == mouse)
                            if len(chosen):
                                ulm_sum[label][mouse] += np.asarray(log_cp[chosen].sum(axis=0)).ravel()
                                ulm_n[label][mouse] += len(chosen)
                        wt_local_rows = np.flatnonzero(discovery)
                        pending_wt_counts[label] = sp.vstack(
                            [pending_wt_counts[label], raw_block[wt_local_rows]], format='csr')
                        pending_wt_meta[label] = pd.concat(
                            [pending_wt_meta[label], frame.iloc[wt_local_rows].reset_index(drop=True)],
                            ignore_index=True)
                        wt_random_rows_written[label] += len(wt_local_rows)
                        while pending_wt_counts[label].shape[0] >= random_chunk_rows:
                            wt_block = pending_wt_counts[label][:random_chunk_rows].tocsr()
                            wt_meta_block = pending_wt_meta[label].iloc[:random_chunk_rows].reset_index(drop=True)
                            _write_mtx_chunk(wt_random_work[label], wt_chunk_number[label],
                                             wt_block, wt_meta_block)
                            wt_chunk_number[label] += 1
                            pending_wt_counts[label] = pending_wt_counts[label][random_chunk_rows:].tocsr()
                            pending_wt_meta[label] = pending_wt_meta[label].iloc[random_chunk_rows:].reset_index(drop=True)
                pending_counts[label] = sp.vstack([pending_counts[label], raw_block], format='csr')
                pending_meta[label] = pd.concat([pending_meta[label], frame], ignore_index=True)
                while pending_counts[label].shape[0] >= chunk_rows:
                    block = pending_counts[label][:chunk_rows].tocsr()
                    block_meta = pending_meta[label].iloc[:chunk_rows].reset_index(drop=True)
                    directory = granule_work if label == GRANULE else target_work[label]
                    _write_mtx_chunk(directory, chunk_number[label], block, block_meta)
                    chunk_number[label] += 1
                    pending_counts[label] = pending_counts[label][chunk_rows:].tocsr()
                    pending_meta[label] = pending_meta[label].iloc[chunk_rows:].reset_index(drop=True)

        for label in (*POPULATIONS, GRANULE):
            directory = granule_work if label == GRANULE else target_work[label]
            if pending_counts[label].shape[0]:
                _write_mtx_chunk(directory, chunk_number[label], pending_counts[label], pending_meta[label].reset_index(drop=True))
                chunk_number[label] += 1
            expected_n = len(metadata_by_label[label])
            assert sum(pd.read_csv(p, sep='\t').shape[0] for p in sorted(directory.glob('chunk_*.nuclei.tsv'))) == expected_n
            metadata_frame = metadata_by_label[label]
            assert np.isfinite(metadata_frame.log10_total_umi.to_numpy(dtype=np.float64)).all()
            assert (metadata_frame.n_detected.to_numpy(dtype=np.int64) >= 0).all()
            metadata_name = (f'{POP_SLUG[label]}_nuclei.tsv' if label in POPULATIONS
                             else 'granule_leak_only_nuclei.tsv')
            metadata_frame.to_csv(work / metadata_name, sep='\t', index=False)
        for population in POPULATIONS:
            if pending_wt_counts[population].shape[0]:
                _write_mtx_chunk(wt_random_work[population], wt_chunk_number[population],
                                 pending_wt_counts[population],
                                 pending_wt_meta[population].reset_index(drop=True))
                wt_chunk_number[population] += 1
            expected_wt_n = int((metadata_by_label[population].genotype.eq('WT')
                                 & metadata_by_label[population].age_weeks.astype(int).isin(DISCOVERY_AGES)).sum())
            observed_wt_n = sum(pd.read_csv(path, sep='\t').shape[0]
                                for path in sorted(wt_random_work[population].glob('chunk_*.nuclei.tsv')))
            assert observed_wt_n == expected_wt_n == wt_random_rows_written[population]
        for p in POPULATIONS:
            assert wt_cell_count[p] > 0 and all(ulm_n[p][m] > 0 for m in mouse_ids[p])
        mean_expression = {p: mean_cp10k_sum[p] / wt_cell_count[p] for p in POPULATIONS}
        ulm_profiles = {p: {m: ulm_sum[p][m] / ulm_n[p][m] for m in mouse_ids[p]} for p in POPULATIONS}
        bergmann_meta = metadata_by_label['Bergmann glia']
        bergmann_wt = bergmann_meta.genotype.eq('WT').to_numpy() & bergmann_meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()
        assert np.isfinite(bergmann_neighbour_score[bergmann_wt]).all()
        pd.DataFrame({'nucleus_id': bergmann_meta.loc[bergmann_wt, 'nucleus_id'].astype(str).to_numpy(),
                      'mouse_id': bergmann_meta.loc[bergmann_wt, 'mouse_id'].astype(str).to_numpy(),
                      'age_weeks': bergmann_meta.loc[bergmann_wt, 'age_weeks'].astype(int).to_numpy(),
                      'neighbour_score': bergmann_neighbour_score[bergmann_wt]}).to_csv(
                          work / 'BERGMANN_T6B_NEIGHBOUR_SCORE.tsv.gz', sep='\t', index=False, compression='gzip')
        details = dict(n_original_gene_columns=int(len(symbols)), n_unique_gene_symbols=int(len(gene_names)),
                       n_duplicate_symbol_columns_aggregated=n_duplicate_symbol_columns_aggregated,
                       n_distinct_duplicated_symbols=int(len(duplicate_symbol_names)),
                       duplicate_symbol_names=duplicate_symbol_names,
                       t6b_neighbour_marker_symbols=sorted(set(neighbor_marker_symbols)),
                       wt_discovery_nuclei_by_population={
                           p: int(wt_random_rows_written[p]) for p in POPULATIONS},
                       wt_discovery_random_chunk_counts={
                           p: int(wt_chunk_number[p]) for p in POPULATIONS},
                       wt_discovery_random_chunk_rows=int(random_chunk_rows),
                       n_t6b_neighbour_markers_present=int(len(marker_indices)))
    return gene_names, mean_expression, ulm_profiles, ulm_n, details


def _validate_env() -> dict[str, Any]:
    explicit = subprocess.run([str(CONDA), 'list', '--explicit', '-p', str(ENV)],
                              check=True, capture_output=True).stdout
    digest = sha256_bytes(explicit)
    assert digest == BASE_ENV_EXPLICIT_SHA256, ('clone explicit list differs from unchanged single_cell base', digest)
    import anndata
    import decoupler as dc
    import h5py
    import numpy as np
    import scipy
    assert sys.version.split()[0] == '3.9.23'
    py_versions = dict(anndata=anndata.__version__, decoupler=dc.__version__, h5py=h5py.__version__,
                       numpy=np.__version__, pandas=pd.__version__, python=sys.version.split()[0],
                       scipy=scipy.__version__)
    expected_py = dict(anndata='0.10.8', decoupler=DECOUPLER_VERSION, h5py='3.13.0',
                       numpy='1.26.4', pandas='2.3.3', python='3.9.23', scipy='1.13.1')
    assert py_versions == expected_py, ('C clone core versions changed', py_versions)
    assert dc.__version__ == DECOUPLER_VERSION
    r_expression = "p<-c('UCell','BiocManager','BiocVersion','BiocParallel','BiocNeighbors','msigdbr'); v<-vapply(p,function(x) as.character(utils::packageVersion(x)),character(1)); cat(paste(c(as.character(getRversion()),v,paste(.libPaths(),collapse=';')),collapse='\\t'))"
    r_out = subprocess.run([str(ENV / 'bin/Rscript'), '--vanilla', '-e', r_expression],
                           check=True, capture_output=True, text=True, env={**os.environ, 'R_LIBS_USER': str(ENV / 'lib/R/library')}).stdout.strip()
    r_fields = r_out.split('\t')
    r_names = ['UCell', 'BiocManager', 'BiocVersion', 'BiocParallel', 'BiocNeighbors', 'msigdbr']
    expected_r = ['2.6.2', '1.30.27', '3.18.1', '1.36.0', '1.20.2', '25.1.1']
    assert len(r_fields) == 8 and r_fields[0] == '4.3.1' and r_fields[1:7] == expected_r
    assert all(str(ENV) in p for p in r_fields[7].split(';')), r_fields[7]
    r_version_file = WORKSPACE / 'env/R_PACKAGE_VERSIONS.tsv'
    assert r_version_file.is_file()
    return dict(clone_path=str(ENV), conda_explicit_sha256=digest,
                conda_explicit_bytes=len(explicit), python_versions=py_versions,
                R=r_fields[0], R_package_versions=dict(zip(r_names, r_fields[1:7])),
                R_package_versions_tsv_sha256=sha256(r_version_file), R_libPaths=r_fields[7],
                conda_explicit_command=f'conda list --explicit -p {ENV}')


def _ucell_invocation(chunk_mtx: Path, genes_path: Path, features_path: Path,
                      metadata_path: Path, random_dir: Optional[Path], out_prefix: Path,
                      include_random: bool, max_ranks: tuple[int, ...] = MAX_RANKS,
                      mode: str = 'standard', base_scores_path: Optional[Path] = None) -> tuple[list[str], dict[str, str]]:
    command = [str(ENV / 'bin/Rscript'), '--vanilla', str(UCELL_HELPER),
               str(chunk_mtx), str(genes_path), str(features_path), str(metadata_path),
               str(random_dir) if include_random and random_dir else '-', str(out_prefix),
               ','.join(map(str, max_ranks)), str(CHUNK_ROWS), mode]
    if mode == 'loo':
        if base_scores_path is None:
            raise ValueError('QC6 LOO mode requires the corresponding observed score block')
        command.append(str(base_scores_path))
    env = dict(os.environ)
    env['R_LIBS_USER'] = str(ENV / 'lib/R/library')
    env['PATH'] = str(ENV / 'bin') + os.pathsep + env.get('PATH', '')
    env['LD_LIBRARY_PATH'] = str(ENV / 'lib') + os.pathsep + env.get('LD_LIBRARY_PATH', '')
    return command, env


def _write_ucell_logs(out_prefix: Path, stdout: str, stderr: str) -> None:
    if stdout.strip():
        with (out_prefix.parent / f'{out_prefix.name}.stdout.txt').open('x', encoding='utf-8') as stream:
            stream.write(stdout)
    # Always retain the step stderr, including an empty file on success.
    with (out_prefix.parent / f'{out_prefix.name}.stderr.txt').open('x', encoding='utf-8') as stream:
        stream.write(stderr)


def _run_ucell_chunk(chunk_mtx: Path, genes_path: Path, features_path: Path,
                     metadata_path: Path, random_dir: Optional[Path], out_prefix: Path,
                     include_random: bool, max_ranks: tuple[int, ...] = MAX_RANKS,
                     mode: str = 'standard', base_scores_path: Optional[Path] = None) -> None:
    command, env = _ucell_invocation(chunk_mtx, genes_path, features_path, metadata_path,
                                      random_dir, out_prefix, include_random, max_ranks,
                                      mode, base_scores_path)
    result = subprocess.run(command, check=True, capture_output=True, text=True, env=env)
    _write_ucell_logs(out_prefix, result.stdout, result.stderr)


def _run_ucell_chunks_parallel(tasks: list[dict[str, Any]],
                               workers: int = DEFAULT_UCELL_PARALLEL_WORKERS) -> None:
    """Run each task once, as local subprocesses, at most `workers` at a time."""
    if workers < 1:
        raise ValueError('UCell parallel worker count must be positive')
    prefixes = [str(task['out_prefix']) for task in tasks]
    if len(prefixes) != len(set(prefixes)):
        raise ValueError('parallel UCell tasks must have unique output prefixes')
    for first in range(0, len(tasks), workers):
        batch = tasks[first:first + workers]
        running = []
        failures = []
        for task in batch:
            out_prefix = Path(task['out_prefix'])
            command, env = _ucell_invocation(**task)
            launch_command = command
            try:
                process = subprocess.Popen(launch_command, stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE, text=True, env=env)
                running.append((process, out_prefix))
            except OSError as exc:
                _write_ucell_logs(out_prefix, '', f'failed to launch task: {exc}\n')
                failures.append((out_prefix, None, str(exc)))
                break
        for process, out_prefix in running:
            stdout, stderr = process.communicate()
            _write_ucell_logs(out_prefix, stdout, stderr)
            if process.returncode:
                failures.append((out_prefix, process.returncode, stderr[-2000:]))
        if failures:
            details = '; '.join(f'{path} exit={code}: {stderr}'
                                for path, code, stderr in failures)
            raise RuntimeError(f'UCell task batch failed without retry: {details}')

def _make_loo_features(records: list[dict], enriched_ids: set[str]) -> tuple[pd.DataFrame, list[str]]:
    by_id = {row['feature_id']: row for row in records}
    base_ids = sorted(enriched_ids)
    rows = []
    for base_i, feature_id in enumerate(base_ids):
        genes = list(by_id[feature_id]['matched_genes'])
        for gene_i, omitted in enumerate(genes):
            loo_id = f'LOO__{feature_id}__{gene_i:03d}'
            retained = genes[:gene_i] + genes[gene_i + 1:]
            rows.append(dict(feature_id=loo_id, base_feature_id=feature_id,
                             base_feature_index=base_i, omitted_gene=omitted,
                             matched_gene_symbols=';'.join(retained)))
    return pd.DataFrame(rows, columns=['feature_id', 'base_feature_id', 'base_feature_index',
                                       'omitted_gene', 'matched_gene_symbols']), base_ids


def qc6_leak_flag(rho: float, driver_genes: set[str], top_neighbour_genes: set[str],
                  cutoff: float = QC6_RHO_CUTOFF) -> bool:
    return bool((np.isfinite(rho) and rho >= cutoff) or bool(driver_genes & top_neighbour_genes))


def run_qc6_bergmann(matrix: np.ndarray, metadata: pd.DataFrame, records: list[dict],
                     enriched_ids: set[str], neighbour_table: pd.DataFrame,
                     top_marker_table: pd.DataFrame, counts_dir: Path,
                     gene_names_path: Path, scratch: Path, out_dir: Path,
                     workers: int = DEFAULT_UCELL_PARALLEL_WORKERS) -> dict[str, Any]:
    """QC6: full-score WT correlation plus per-gene leave-one-out driver attribution."""
    pop = 'Bergmann glia'
    slug = POP_SLUG[pop]
    enriched_ids = set(enriched_ids)
    top_markers = set(top_marker_table.loc[
        top_marker_table.population.isin([GRANULE, 'Purkinje cell']), 'symbol'].astype(str))
    neighbour = neighbour_table.set_index('nucleus_id').neighbour_score
    metadata_ids = metadata.nucleus_id.astype(str)
    wt_rows = np.flatnonzero(metadata.genotype.astype(str).str.upper().eq('WT').to_numpy()
                             & metadata.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy())
    wt_nuclei = metadata_ids.iloc[wt_rows].tolist()
    neighbour_scores = neighbour.reindex(wt_nuclei).to_numpy(dtype=np.float64)
    if not np.isfinite(neighbour_scores).all():
        raise ValueError('QC6 T6b neighbour score is not aligned/finite for WT Bergmann nuclei')
    loo_table, base_ids = _make_loo_features(records, enriched_ids)
    if not base_ids:
        empty_flags = pd.DataFrame(columns=['population', 'feature_id', 'set_name', 'WT_Bergmann_spearman_rho_to_T6b_neighbour_score',
                                             'n_set_genes_overlapping_top_Granule_or_Purkinje', 'n_LOO_driver_genes_overlapping_top_markers',
                                             'QC6_driver_marker_warning', 'QC6_neighbour_correlation_warning',
                                             'QC6_warning_conditions', 'flagged_by_QC6'])
        empty_flags.to_csv(out_dir / 'QC6_BERGMANN_LEAK_FLAGS.tsv', sep='\t', index=False)
        pd.DataFrame(columns=['feature_id', 'base_feature_id', 'omitted_gene', 'mean_absolute_LOO_score_change',
                              'driver_rank', 'driver_top_decile', 'top_Granule_or_Purkinje_marker']).to_csv(
                                  out_dir / 'QC6_BERGMANN_LOO_DRIVER_DETAILS.tsv.gz', sep='\t', index=False, compression='gzip')
        return dict(n_primary_enriched_sets=0, n_driver_rows=0, n_flagged_sets=0,
                    T6b_script_sha256=T6B_SCRIPT_SHA256, flagged_feature_ids=[],
                    warning_conditions_by_feature={})

    loo_abs_sum = np.zeros(len(loo_table), dtype=np.float64)
    loo_wt_n = 0
    cursor = 0
    all_feature_ids = [row['feature_id'] for row in records]
    score_col = {feature: i for i, feature in enumerate(all_feature_ids)}
    chunk_mtx_files = sorted(counts_dir.glob('chunk_*.mtx'))
    # Cache chunk metadata and offsets once, then batch enriched sets so retained
    # UCell output remains bounded even when many sets enrich.
    chunk_specs = []
    cursor = 0
    for chunk_mtx in chunk_mtx_files:
        prefix = chunk_mtx.with_suffix('')
        meta_path = Path(str(prefix) + '.nuclei.tsv')
        chunk_meta = pd.read_csv(meta_path, sep='\t', dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str})
        n_chunk = len(chunk_meta)
        local_wt = (chunk_meta.genotype.astype(str).str.upper().eq('WT').to_numpy()
                    & chunk_meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy())
        chunk_specs.append((chunk_mtx, meta_path, cursor, n_chunk, local_wt))
        cursor += n_chunk
    if cursor != len(metadata):
        raise ValueError('QC6 chunk traversal differs from the Bergmann score-matrix row universe')
    for first_base in range(0, len(base_ids), QC6_LOO_BASE_SET_BATCH):
        batch_base_ids = base_ids[first_base:first_base + QC6_LOO_BASE_SET_BATCH]
        batch_base_index = {feature_id: i for i, feature_id in enumerate(batch_base_ids)}
        batch_rows = np.flatnonzero(loo_table.base_feature_id.astype(str).isin(batch_base_ids).to_numpy())
        batch_features = loo_table.iloc[batch_rows].copy()
        batch_features['base_feature_index'] = batch_features.base_feature_id.astype(str).map(batch_base_index).astype(int)
        loo_path = scratch / f'QC6_BERGMANN_LOO_FEATURES_{first_base:06d}.tsv'
        batch_features.to_csv(loo_path, sep='\t', index=False)
        base_cols = [score_col[feature] for feature in batch_base_ids]
        loo_tasks = []
        loo_specs = []
        for chunk_mtx, meta_path, row_start, n_chunk, local_wt in chunk_specs:
            if not local_wt.any():
                continue
            base_block = np.asarray(matrix[row_start:row_start + n_chunk, base_cols], dtype=np.float64)
            base_file = scratch / f'qc6_base_scores_{first_base:06d}_{row_start:08d}.bin'
            base_block.astype('<f8', copy=False).ravel(order='F').tofile(base_file)
            loo_prefix = scratch / f'qc6_loo_{first_base:06d}_{row_start:08d}'
            loo_tasks.append(dict(
                chunk_mtx=chunk_mtx, genes_path=gene_names_path, features_path=loo_path,
                metadata_path=meta_path, random_dir=None, out_prefix=loo_prefix,
                include_random=False, max_ranks=(PRIMARY_MAX_RANK,), mode='loo',
                base_scores_path=base_file))
            loo_specs.append((loo_prefix, int(local_wt.sum())))
        _run_ucell_chunks_parallel(loo_tasks, workers)
        for loo_prefix, n_local_wt in loo_specs:
            effect_path = Path(f'{loo_prefix}.loo_abs_delta_sum.bin')
            effect = np.fromfile(effect_path, dtype='<f8')
            if effect.size != len(batch_rows):
                raise ValueError('QC6 LOO effect vector has the wrong length')
            loo_abs_sum[batch_rows] += effect
            if first_base == 0:
                loo_wt_n += n_local_wt
        if cursor != len(metadata) or loo_wt_n != len(wt_rows):
            raise ValueError('QC6 chunk traversal differs from the WT Bergmann row universe')
    loo_table['mean_absolute_LOO_score_change'] = loo_abs_sum / loo_wt_n

    record_by_id = {row['feature_id']: row for row in records}
    detail_frames = []
    flag_rows = []
    for feature_id in base_ids:
        row = record_by_id[feature_id]
        j = score_col[feature_id]
        full_score = np.asarray(matrix[wt_rows, j], dtype=np.float64)
        rho = float(spearmanr(full_score, neighbour_scores).statistic)
        subset = loo_table.loc[loo_table.base_feature_id.eq(feature_id)].copy()
        subset = subset.sort_values(['mean_absolute_LOO_score_change', 'omitted_gene'],
                                    ascending=[False, True], kind='mergesort').reset_index(drop=True)
        subset['driver_rank'] = np.arange(1, len(subset) + 1)
        n_driver = max(1, int(np.ceil(len(subset) * QC6_LOO_DRIVER_TOP_FRACTION)))
        subset['driver_top_decile'] = subset.driver_rank <= n_driver
        subset['top_Granule_or_Purkinje_marker'] = subset.omitted_gene.astype(str).isin(top_markers)
        subset['WT_Bergmann_spearman_rho_to_T6b_neighbour_score'] = rho
        subset['primary_set_score_vs_neighbour_rho_pass'] = bool(np.isfinite(rho) and rho >= QC6_RHO_CUTOFF)
        driver_genes = set(subset.loc[subset.driver_top_decile, 'omitted_gene'].astype(str))
        set_markers = set(row['matched_genes']) & top_markers
        driver_markers = driver_genes & top_markers
        driver_marker_warning = bool(driver_markers)
        neighbour_correlation_warning = bool(np.isfinite(rho) and rho >= QC6_RHO_CUTOFF)
        flagged = qc6_leak_flag(rho, driver_genes, top_markers)
        warning_conditions = []
        if driver_marker_warning:
            warning_conditions.append('LOO_driver_overlaps_T6b_Granule_or_Purkinje_top_markers')
        if neighbour_correlation_warning:
            warning_conditions.append('WT_Bergmann_set_score_rho_ge_0.4_to_neighbour_score')
        flag_rows.append(dict(population=pop, feature_id=feature_id, set_name=row['set_name'],
            n_matched_genes=len(row['matched_genes']),
            WT_Bergmann_spearman_rho_to_T6b_neighbour_score=rho,
            n_set_genes_overlapping_top_Granule_or_Purkinje=len(set_markers),
            set_marker_overlap_genes=';'.join(sorted(set_markers)),
            n_LOO_driver_genes=len(driver_genes), n_LOO_driver_genes_overlapping_top_markers=len(driver_markers),
            LOO_driver_marker_genes=';'.join(sorted(driver_markers)),
            rho_cutoff=QC6_RHO_CUTOFF, driver_definition='top 10% by mean absolute UCell-score change after leave-one-gene-out',
            QC6_driver_marker_warning=driver_marker_warning,
            QC6_neighbour_correlation_warning=neighbour_correlation_warning,
            QC6_warning_conditions=';'.join(warning_conditions),
            flagged_by_QC6=flagged))
        subset.insert(0, 'set_name', row['set_name'])
        detail_frames.append(subset)
    flags = pd.DataFrame(flag_rows)
    flags.to_csv(out_dir / 'QC6_BERGMANN_LEAK_FLAGS.tsv', sep='\t', index=False)
    pd.concat(detail_frames, ignore_index=True).to_csv(
        out_dir / 'QC6_BERGMANN_LOO_DRIVER_DETAILS.tsv.gz', sep='\t', index=False, compression='gzip')
    return dict(n_primary_enriched_sets=len(base_ids), n_driver_rows=int(len(loo_table)),
                n_flagged_sets=int(flags.flagged_by_QC6.sum()), flagged_feature_ids=flags.loc[
                    flags.flagged_by_QC6, 'feature_id'].astype(str).tolist(),
                warning_conditions_by_feature={str(row.feature_id): str(row.QC6_warning_conditions).split(';')
                                               if str(row.QC6_warning_conditions) else []
                                               for row in flags.itertuples()},
                T6b_script_sha256=T6B_SCRIPT_SHA256, WT_Bergmann_nuclei=len(wt_rows),
                rho_cutoff=QC6_RHO_CUTOFF, LOO_base_set_batch_size=QC6_LOO_BASE_SET_BATCH,
                driver_definition='top 10% per enriched set by mean absolute UCell-score change after leave-one-gene-out; compared with T6b top-50 Granule/Purkinje WT markers')


def _assert_wt_random_chunk_nucleus_order(expected_nucleus_ids: Any,
                                          chunk_nucleus_ids: Any) -> None:
    """Require WT-only random chunks to preserve the full metadata row order."""
    expected = [str(value) for value in expected_nucleus_ids]
    observed = [str(value) for value in chunk_nucleus_ids]
    assert observed == expected, (
        'WT random-chunk nucleus_id order differs from metadata.nucleus_id.iloc[wt_row_ids]')


def _compute_ulm(ulm_profiles: dict[str, np.ndarray], mouse_ids: list[str],
                 records: list[dict], gene_names: np.ndarray) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    import decoupler as dc
    genes = list(map(str, gene_names))
    feature_ids = [str(row['feature_id']) for row in records]
    mat = pd.DataFrame(index=mouse_ids, columns=genes, dtype=np.float64)
    gene_idx = {g: i for i, g in enumerate(genes)}
    # Each mouse profile is already the mean of nucleus-level log1p(CP10k) values.
    for i, mouse in enumerate(mouse_ids):
        profile = ulm_profiles[mouse]
        mat.loc[mouse, genes] = [profile[gene_idx[g]] for g in genes]
    net_rows = []
    for row in records:
        net_rows.extend(dict(source=row['feature_id'], target=g, weight=1.0) for g in row['matched_genes'])
    network = pd.DataFrame(net_rows, columns=['source', 'target', 'weight'])
    estimate, _pvals = dc.run_ulm(mat=mat, net=network, source='source', target='target', weight='weight',
                                  min_n=MIN_MATCHED_GENES, verbose=False, use_raw=False)
    returned_ids = list(map(str, estimate.columns))
    if len(returned_ids) != len(set(returned_ids)):
        raise ValueError('decoupler ULM returned duplicate source IDs')
    unexpected_ids = sorted(set(returned_ids) - set(feature_ids))
    if unexpected_ids:
        raise ValueError(f'decoupler ULM returned unknown source IDs: {unexpected_ids}')
    returned_id_set = set(returned_ids)
    missing_ids = [feature_id for feature_id in feature_ids if feature_id not in returned_id_set]
    estimate = estimate.copy()
    estimate.columns = returned_ids
    estimate = estimate.reindex(index=mouse_ids, columns=feature_ids)
    # decoupler drops sources below min_n after removing all-zero assay genes. A wholly
    # absent source is ULM-missing, never ULM-positive; encode it as zero before enrichment.
    if missing_ids:
        estimate.loc[:, missing_ids] = 0.0
    values = estimate.to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('decoupler ULM returned non-finite estimates for returned sources or mice')
    return estimate, values, missing_ids


def _mouse_ids(meta: pd.DataFrame) -> list[str]:
    wt = meta.genotype.eq('WT') & meta.age_weeks.isin(DISCOVERY_AGES)
    result = sorted(meta.loc[wt, 'mouse_id'].astype(str).unique())
    if len(result) != EXPECTED_MICE_PER_POP:
        raise ValueError(f'expected 12 WT discovery mice, found {len(result)}: {result}')
    return result


def _summarize_scores(population: str, matrix: np.ndarray, meta: pd.DataFrame,
                      records: list[dict], ulm_values: np.ndarray,
                      mouse_median_shares_by_rank: dict[int, np.ndarray], out_dir: Path) -> dict[str, Any]:
    feature_ids = [r['feature_id'] for r in records]
    n_features = len(records)
    mouse_ids = _mouse_ids(meta)
    ulm_frame = pd.DataFrame(ulm_values, index=mouse_ids, columns=feature_ids)
    ulm_frame.index.name = 'mouse_id'
    ulm_frame.to_csv(out_dir / f'{POP_SLUG[population]}_ULM_BY_MOUSE.tsv', sep='\t')

    enrichment_frame, enriched = summarize_enrichment(mouse_median_shares_by_rank, ulm_values, mouse_ids, feature_ids)
    enrichment_frame.insert(0, 'population', population)
    enrichment_frame.to_csv(out_dir / f'{POP_SLUG[population]}_ENRICHMENT_SHARES.tsv.gz', sep='\t', index=False, compression='gzip')
    enriched_rows = []
    for rank in MAX_RANKS:
        shares = mouse_median_shares_by_rank[rank]
        for feature_id in sorted(enriched[rank]):
            j = feature_ids.index(feature_id)
            enriched_rows.append(dict(population=population, maxRank=rank,
                feature_id=feature_id, set_name=records[j]['set_name'],
                n_matched_genes=records[j]['n_matched_genes'],
                n_WT_mice_median_share_ge_0_95=int((shares[:, j] >= ENRICHMENT_SHARE_CUTOFF).sum()),
                n_WT_mice_median_share_and_ULM_positive=int(((shares[:, j] >= ENRICHMENT_SHARE_CUTOFF) & (ulm_values[:, j] > 0)).sum())))
    pd.DataFrame(enriched_rows).to_csv(out_dir / f'{POP_SLUG[population]}_ENRICHED_SETS.tsv', sep='\t', index=False)

    ceiling_share = np.zeros(n_features, dtype=np.int64)
    near_ceiling_share = np.zeros(n_features, dtype=np.int64)
    zero_share = np.zeros(n_features, dtype=np.int64)
    ge099_per_nucleus = np.empty(len(matrix), dtype=np.int64)
    exact_ceiling_per_nucleus = np.empty(len(matrix), dtype=np.int64)
    iqr_per_nucleus = np.empty(len(matrix), dtype=np.float64)
    for first in range(0, len(matrix), CHUNK_ROWS):
        last = min(first + CHUNK_ROWS, len(matrix))
        block = np.asarray(matrix[first:last], dtype=np.float64)
        ceiling_share += np.sum(block >= 1.0 - 1e-12, axis=0)
        near_ceiling_share += np.sum(block >= 0.99, axis=0)
        zero_share += np.sum(block <= 1e-12, axis=0)
        ge099_per_nucleus[first:last] = np.sum(block >= 0.99, axis=1)
        exact_ceiling_per_nucleus[first:last] = np.sum(block >= 1.0 - 1e-12, axis=1)
        quartiles = np.quantile(block, [0.25, 0.75], axis=1)
        iqr_per_nucleus[first:last] = quartiles[1] - quartiles[0]
    ceiling_share = ceiling_share / len(matrix)
    near_ceiling_share = near_ceiling_share / len(matrix)
    zero_share = zero_share / len(matrix)
    iqr_frame = pd.DataFrame({'population': population, 'nucleus_id': meta.nucleus_id.astype(str),
                              'within_nucleus_score_IQR': iqr_per_nucleus,
                              'n_sets_score_ge_0_99': ge099_per_nucleus,
                              'n_sets_score_exactly_1': exact_ceiling_per_nucleus})
    iqr_frame.to_csv(out_dir / f'{POP_SLUG[population]}_WITHIN_NUCLEUS_IQR.tsv.gz', sep='\t', index=False, compression='gzip')

    gate_summary = qc_hard_gates(matrix, float(np.mean(iqr_per_nucleus == 0.0)), ge099_per_nucleus)
    covariates = {'log10_total_umi': meta.log10_total_umi.to_numpy(dtype=float),
                  'n_detected': meta.n_detected.to_numpy(dtype=float)}
    corr_rows = []
    for j, row in enumerate(records):
        score = np.asarray(matrix[:, j], dtype=np.float64)
        for covariate, values in covariates.items():
            rho = float(spearmanr(score, values).statistic)
            corr_rows.append(dict(population=population, feature_id=row['feature_id'], set_name=row['set_name'],
                                  covariate=covariate, spearman_rho=rho,
                                  abs_rho_flag=bool(np.isfinite(rho) and abs(rho) > CORRELATION_FLAG)))
    corr_frame = pd.DataFrame(corr_rows)
    corr_frame.to_csv(out_dir / f'{POP_SLUG[population]}_SCORE_COVARIATE_CORRELATIONS.tsv', sep='\t', index=False)

    full_means, down_means = downsample_mouse_means(matrix, meta, mouse_ids, seed=SEED)
    r_rows = []
    for j, row in enumerate(records):
        if np.std(full_means[:, j]) == 0 or np.std(down_means[:, j]) == 0:
            r = np.nan
        else:
            r = float(np.corrcoef(full_means[:, j], down_means[:, j])[0, 1])
        r_rows.append(dict(population=population, feature_id=row['feature_id'], set_name=row['set_name'],
                           seed=SEED, sampled_fraction=0.5, Pearson_r_full_vs_downsampled_mouse_means=r,
                           flag_r_below_0_8=bool(not np.isfinite(r) or r < QC4_DOWNSAMPLE_R_FLAG)))
    pd.DataFrame(r_rows).to_csv(out_dir / f'{POP_SLUG[population]}_WT_DOWNSAMPLE_50PCT_R.tsv', sep='\t', index=False)

    wt_mask = meta.genotype.astype(str).str.upper().eq('WT').to_numpy() & meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()
    primary_ids = sorted(enriched[PRIMARY_MAX_RANK])
    if primary_ids:
        primary_cols = [feature_ids.index(feature_id) for feature_id in primary_ids]
        wt_rows = np.flatnonzero(wt_mask)
        wt_score_slice = np.asarray(matrix[np.ix_(wt_rows, primary_cols)], dtype=np.float64)
    else:
        wt_score_slice = np.empty((int(wt_mask.sum()), 0), dtype=np.float64)
    primary_groups, primary_reps = group_redundant_sets(enriched[PRIMARY_MAX_RANK], records,
                                                         wt_score_slice, primary_ids)
    primary_groups.insert(0, 'population', population)
    primary_reps.insert(0, 'population', population)
    primary_groups.to_csv(out_dir / f'{POP_SLUG[population]}_ENRICHED_GROUPS.tsv', sep='\t', index=False)
    primary_reps.to_csv(out_dir / f'{POP_SLUG[population]}_GROUP_REPRESENTATIVES.tsv', sep='\t', index=False)

    jac_rows = []
    for a, b in ((1000, 1500), (1500, 3000), (1000, 3000)):
        sa, sb = enriched[a], enriched[b]
        jac = len(sa & sb) / len(sa | sb) if sa | sb else 1.0
        jac_rows.append(dict(population=population, maxRank_a=a, maxRank_b=b,
                             enriched_n_a=len(sa), enriched_n_b=len(sb),
                             intersection_n=len(sa & sb), union_n=len(sa | sb), jaccard=jac,
                             flag_jaccard_below_0_7=qc5_jaccard_flag(jac)))
    pd.DataFrame(jac_rows).to_csv(out_dir / f'{POP_SLUG[population]}_MAXRANK_ENRICHED_LIST_JACCARD.tsv',
                                  sep='\t', index=False)

    qc_flags = []
    if corr_frame.abs_rho_flag.any():
        qc_flags.append('QC3: score correlation |rho|>0.4 with log10_total_umi or n_detected')
    if any(row['flag_r_below_0_8'] for row in r_rows):
        qc_flags.append('QC4: one or more fixed-seed 50% WT downsample correlations are below 0.8 or undefined')
    if any(row['flag_jaccard_below_0_7'] for row in jac_rows):
        qc_flags.append('QC5: one or more maxRank enriched-list Jaccards are below 0.7')
    pd.DataFrame([dict(population=population, feature_id=records[j]['feature_id'], set_name=records[j]['set_name'],
                       exact_ceiling_score_share=float(ceiling_share[j]), near_ceiling_score_share=float(near_ceiling_share[j]),
                       zero_score_share=float(zero_share[j]),
                       any_ceiling=bool(ceiling_share[j] > 0), all_ceiling=bool(ceiling_share[j] >= 1),
                       any_near_ceiling=bool(near_ceiling_share[j] > 0), all_near_ceiling=bool(near_ceiling_share[j] >= 1),
                       any_zero=bool(zero_share[j] > 0), all_zero=bool(zero_share[j] >= 1))
                  for j in range(n_features)]).to_csv(out_dir / f'{POP_SLUG[population]}_CEILING_ZERO_SHARES.tsv', sep='\t', index=False)
    qc_summary = dict(population=population, n_features=n_features,
        sets_with_any_exact_ceiling_score=int(np.sum(ceiling_share > 0)),
        sets_fully_at_exact_ceiling=int(np.sum(ceiling_share >= 1.0)),
        sets_with_any_near_ceiling_score=int(np.sum(near_ceiling_share > 0)),
        sets_fully_at_or_above_0_99=int(np.sum(near_ceiling_share >= 1.0)),
        sets_with_any_zero_score=int(np.sum(zero_share > 0)), sets_fully_zero=int(np.sum(zero_share >= 1.0)),
        maximum_exact_ceiling_share=float(np.max(ceiling_share)),
        maximum_near_ceiling_share=float(np.max(near_ceiling_share)), maximum_zero_share=float(np.max(zero_share)),
        median_within_nucleus_IQR=float(np.median(iqr_per_nucleus)),
        score_covariate_flag_count=int(corr_frame.abs_rho_flag.sum()),
        qc3_to_qc6_flags=qc_flags, hard_gates=gate_summary,
        enriched_count_by_maxRank={str(k): len(v) for k, v in enriched.items()})
    return dict(qc_summary=qc_summary, enriched=enriched, rank_jaccards=jac_rows,
                mouse_ids=mouse_ids, full_mouse_means=full_means,
                feature_ids=feature_ids, wt_mask=wt_mask)


def run_real(args: argparse.Namespace) -> dict[str, Any]:
    # The authorization check is the first real-run action; no source file is opened before it.
    if sys.flags.isolated != 1 or not sys.dont_write_bytecode:
        raise RuntimeError('refusing real-data entry point unless launched with python -I -B')
    gate = check_go(args.go_file)
    if args.out.exists():
        raise FileExistsError(f'refusing to overwrite an output or retry a failed run: {args.out}')
    args.out.mkdir(parents=True)
    run_started = time.time()
    try:
        env_record = _validate_env()
        gmt_digest = sha256(args.gmt)
        assert gmt_digest == GMT_SHA256
        for label, path in OFFICIAL_GMTS.items():
            assert path.is_file(), f'missing official GMT for comparison: {label} {path}'
        comparison = compare_official_gmts(args.gmt, OFFICIAL_GMTS)
        if not comparison['exact_names_and_members_equal']:
            _write_json(args.out / 'MSIGDB_COMPARISON.json', comparison)
            raise RuntimeError('native msigdbr GMT differs from official mouse GMT; 712fd1c0 remains pinned, review differences before any analysis')
        pins = {str(path): sha256(path) for path in (
            args.h5ad, args.manifest, args.animal_manifest, args.scoping, args.gmt)}
        assert pins[str(args.h5ad)] == H5AD_SHA256
        assert pins[str(args.manifest)] == NUCLEUS_MANIFEST_SHA256
        assert pins[str(args.animal_manifest)] == ANIMAL_MANIFEST_SHA256
        assert pins[str(args.scoping)] == SCOPING_SHA256
        assert pins[str(args.gmt)] == GMT_SHA256
        nucleus_manifest = pd.read_csv(args.manifest, sep='\t', dtype={'animal_id': str}, low_memory=False)
        animal_manifest = pd.read_csv(args.animal_manifest, sep='\t', dtype={'animal_id': str}, low_memory=False)
        scoping_record = json.loads(args.scoping.read_text(encoding='utf-8'))
        p_meta = build_p_nucleus_metadata(args.h5ad, nucleus_manifest, animal_manifest, scoping_record)
        t6b_wt_cells = _read_t6b_wt_manifest(args.manifest, args.animal_manifest)
        t6b_marker_table, t6b_marker_symbols = compute_t6b_top_markers(args.h5ad, t6b_wt_cells)
        granule_marker_count = int((t6b_marker_table.population == GRANULE).sum())
        granule_marker_count_at_ratio4 = int(((t6b_marker_table.population == GRANULE)
                                              & (t6b_marker_table.own_vs_max_other_mean_ratio >= T6B_PRIMARY_RATIO)).sum())
        purkinje_marker_count = int((t6b_marker_table.population == 'Purkinje cell').sum())
        bergmann_marker_count = int((t6b_marker_table.population == 'Bergmann glia').sum())
        if granule_marker_count != 6:
            raise RuntimeError(f'predeclared QC6 T6b Granule ratio-3 marker count was expected to be 6; found {granule_marker_count}')
        if granule_marker_count_at_ratio4 != 0:
            raise RuntimeError(f'predeclared QC6 T6b Granule ratio-4 marker count was expected to be 0; found {granule_marker_count_at_ratio4}')
        bergmann_self_markers = set(t6b_marker_symbols['Bergmann glia'])
        neighbour_symbols = t6b_bergmann_neighbour_symbols(t6b_marker_symbols)
        if not neighbour_symbols:
            raise RuntimeError('T6b Granule/Purkinje marker union is empty after Bergmann self-marker removal')
        granule = t6b_wt_cells.loc[t6b_wt_cells.cell_type_raw.eq(GRANULE),
                                   ['nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex',
                                     'log10_total_umi', 'n_detected']].copy().reset_index(drop=True)
        scratch = args.out / '.scratch'
        scratch.mkdir()
        gene_names, mean_expression, ulm_profiles, _ulm_n, source_details = stage_source(
            args.h5ad, p_meta, granule, neighbour_symbols, scratch / 'staged_counts',
            random_chunk_rows=args.random_chunk_rows)
        t6b_marker_table.to_csv(args.out / 'T6B_WT_TOP_MARKERS.tsv', sep='\t', index=False)
        _write_json(args.out / 'T6B_MARKER_DEFINITION.json', dict(
            source_script=str(T6B_SCRIPT), source_script_sha256=T6B_SCRIPT_SHA256,
            source_run_folder_read=False, active_unit_floor=T6B_ACTIVE_UNIT_FLOOR,
            WT_ages=list(DISCOVERY_AGES), mean_CP10k_floor=T6B_MIN_MEAN_CP10K,
            detection_prevalence_floor=T6B_MIN_PREVALENCE,
            ratio_thresholds_by_population={GRANULE: T6B_GRANULE_QC6_RATIO,
                                            'Purkinje cell': T6B_PRIMARY_RATIO,
                                            'Bergmann glia': T6B_PRIMARY_RATIO},
            selected_marker_counts_by_population={GRANULE: granule_marker_count,
                                                  'Purkinje cell': purkinje_marker_count,
                                                  'Bergmann glia': bergmann_marker_count},
            granule_markers_at_ratio_4=granule_marker_count_at_ratio4,
            granule_markers_at_ratio_2=49,
            granule_threshold_reason='QC6-only predeclared relaxation: T6b showed 0 Granule markers at ratio 4, 6 at ratio 3, and 49 at ratio 2; ratio 3 was selected as the first relaxation yielding the predeclared six markers.',
            top_n_per_population=T6B_TOP_N,
            marker_score='per-nucleus sum of log1p(CP10k) over the unique-symbol union of top Granule (ratio 3 for QC6 only) and Purkinje (ratio 4) markers after removing Bergmann self-marker symbols, matching T6b self-marker exclusion',
            selected_granule_symbols=t6b_marker_symbols[GRANULE],
            selected_purkinje_symbols=t6b_marker_symbols['Purkinje cell'],
            excluded_Bergmann_self_marker_symbols=sorted(bergmann_self_markers),
            final_neighbour_score_symbols=neighbour_symbols))
        (scratch / 'gene_names.txt').write_text('\n'.join(gene_names) + '\n', encoding='utf-8')
        gmt_sets = parse_gmt(args.gmt)
        records, features = build_matched_sets(gmt_sets, gene_names)
        scoring_plan = {}
        for population in POPULATIONS:
            n_all = len(p_meta[population])
            n_wt = source_details['wt_discovery_nuclei_by_population'][population]
            scoring_plan[population] = dict(
                matched_set_count=len(records),
                all_nuclei_count=n_all,
                all_nuclei_chunk_count=(n_all + CHUNK_ROWS - 1) // CHUNK_ROWS,
                WT_discovery_nuclei_count=n_wt,
                WT_discovery_chunk_count=source_details['wt_discovery_random_chunk_counts'][population],
                WT_random_chunk_rows=args.random_chunk_rows,
                random_feature_batch_count=(len(records) + RANDOM_BATCH_FEATURES - 1) // RANDOM_BATCH_FEATURES,
                maxRanks=list(MAX_RANKS), random_controls_per_set=N_RANDOM_SETS)
            scoring_plan[population]['expected_retained_UCell_binary_bytes'] = int(
                (n_all + n_wt) * len(records) * 8
                + 2 * n_wt * len(records) * len(MAX_RANKS))
        scoring_plan['execution'] = dict(
            observed_scores='all nuclei, maxRank 1500 only, scored once in 256-row population chunks',
            random_controls='WT discovery nuclei only; random scoring is split over nucleus chunks and 128-set control batches',
            random_chunk_rows=args.random_chunk_rows,
            random_feature_batch_size=RANDOM_BATCH_FEATURES,
            parallel_workers=args.workers,
            backend='local subprocesses, at most --workers at a time',
            Rscript_processes='one independent Rscript --vanilla process per chunk/control-batch task; UCell ncores=1')
        population_qc = {}
        population_analysis = {}
        ulm_missing_by_population = {}
        population_features = {}
        for pop in POPULATIONS:
            bins = expression_bin_ids(mean_expression[pop], gene_names)
            pop_slug = POP_SLUG[pop]
            features_path = args.out / f'state_P_{pop_slug}_features.tsv'
            features.to_csv(features_path, sep='\t', index=False)
            population_features[pop] = features
            pop_out = scratch / pop_slug
            pop_out.mkdir()
            random_dir = pop_out / 'random_gmt'
            write_random_batches_for_population(records, gene_names, bins, pop, random_dir)
            meta = p_meta[pop].copy()
            # P independently reconstructs A-compatible metadata; shared readout aligns by nucleus_id.
            obs_meta_path = scratch / 'staged_counts' / f'{pop_slug}_nuclei.tsv'
            metadata = pd.read_csv(obs_meta_path, sep='\t', dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str, 'sex': str})
            metadata['age_weeks'] = metadata.age_weeks.astype(int)
            metadata['genotype'] = metadata.genotype.astype(str).str.upper()
            metadata['mouse_id'] = metadata.mouse_id.astype(str)
            metadata.to_csv(args.out / f'state_P_{pop_slug}_nuclei.tsv', sep='\t', index=False)
            mouse_ids = _mouse_ids(metadata)
            ulm_estimate, ulm_values, ulm_missing_ids = _compute_ulm(
                ulm_profiles[pop], mouse_ids, records, gene_names)
            ulm_missing_by_population[pop] = dict(
                count=len(ulm_missing_ids), feature_ids=ulm_missing_ids)
            score_matrix = np.lib.format.open_memmap(
                args.out / f'state_P_{pop_slug}.npy', mode='w+', dtype=np.float64,
                shape=(len(metadata), len(records)))
            wt_discovery = (metadata.genotype.astype(str).str.upper().eq('WT').to_numpy()
                            & metadata.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy())
            wt_row_ids = np.flatnonzero(wt_discovery)
            share_memmaps = {}
            share_paths = {}
            for rank in MAX_RANKS:
                share_paths[rank] = pop_out / f'per_nucleus_sharecount_{rank}.u8'
                share_memmaps[rank] = np.memmap(share_paths[rank], mode='w+', dtype=np.uint8,
                                                shape=(len(wt_row_ids), len(records)))
            chunks = sorted((scratch / 'staged_counts' / pop_slug).glob('chunk_*.mtx'))
            cursor = 0
            observed_specs = []
            observed_tasks = []
            for chunk_mtx in chunks:
                prefix = chunk_mtx.with_suffix('')
                chunk_meta_path = Path(str(prefix) + '.nuclei.tsv')
                chunk_meta = pd.read_csv(chunk_meta_path, sep='\t', dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str, 'sex': str})
                n_chunk = len(chunk_meta)
                score_prefix = pop_out / f'observed_{cursor:08d}'
                observed_tasks.append(dict(
                    chunk_mtx=chunk_mtx, genes_path=scratch / 'gene_names.txt',
                    features_path=features_path, metadata_path=chunk_meta_path,
                    random_dir=None, out_prefix=score_prefix, include_random=False,
                    max_ranks=(PRIMARY_MAX_RANK,)))
                observed_specs.append((score_prefix, cursor, n_chunk))
                cursor += n_chunk
            if cursor != len(metadata):
                raise ValueError(f'{pop}: observed chunk traversal differs from the P nucleus table')
            _run_ucell_chunks_parallel(observed_tasks, args.workers)
            for score_prefix, row_start, n_chunk in observed_specs:
                primary_path = Path(f'{score_prefix}.score_{PRIMARY_MAX_RANK}.bin')
                score = np.fromfile(primary_path, dtype='<f8')
                assert score.size == n_chunk * len(records)
                score = score.reshape((n_chunk, len(records)), order='F')
                assert np.isfinite(score).all()
                score_matrix[row_start:row_start + n_chunk] = score

            wt_random_dir = scratch / 'staged_counts' / f'{pop_slug}_wt_discovery_random'
            wt_random_chunks = sorted(wt_random_dir.glob('chunk_*.mtx'))
            random_gmt_batches = sorted(random_dir.glob('batch_*.gmt'))
            if len(random_gmt_batches) != (len(records) + RANDOM_BATCH_FEATURES - 1) // RANDOM_BATCH_FEATURES:
                raise ValueError(f'{pop}: random-control GMT batch count differs from feature batch plan')
            random_features_dir = pop_out / 'random_features'
            random_features_dir.mkdir()
            random_feature_batches = []
            for batch_i, gmt_path in enumerate(random_gmt_batches):
                feature_start = batch_i * RANDOM_BATCH_FEATURES
                batch_features = features.iloc[feature_start:feature_start + RANDOM_BATCH_FEATURES].copy()
                if batch_features.empty:
                    raise ValueError(f'{pop}: empty random feature batch {batch_i}')
                batch_features_path = random_features_dir / f'batch_{batch_i:05d}.tsv'
                batch_features.to_csv(batch_features_path, sep='\t', index=False)
                random_feature_batches.append((gmt_path, batch_features_path, feature_start, len(batch_features)))
            random_tasks = []
            random_specs = []
            wt_random_nucleus_ids = []
            for random_chunk_i, chunk_mtx in enumerate(wt_random_chunks):
                prefix = chunk_mtx.with_suffix('')
                chunk_meta_path = Path(str(prefix) + '.nuclei.tsv')
                chunk_meta = pd.read_csv(chunk_meta_path, sep='\t', dtype={
                    'nucleus_id': str, 'mouse_id': str, 'genotype': str, 'sex': str})
                n_chunk = len(chunk_meta)
                wt_random_nucleus_ids.extend(chunk_meta.nucleus_id.astype(str).tolist())
                if not (chunk_meta.genotype.astype(str).str.upper().eq('WT')
                        & chunk_meta.age_weeks.astype(int).isin(DISCOVERY_AGES)).all():
                    raise ValueError(f'{pop}: WT random chunk contains a non-WT discovery nucleus')
                score_prefix = pop_out / f'wt_random_{random_chunk_i:08d}'
                batch_outputs = []
                for batch_i, (batch_gmt, batch_features_path, feature_start, feature_count) in enumerate(random_feature_batches):
                    batch_prefix = pop_out / f'wt_random_{random_chunk_i:05d}_set_batch_{batch_i:05d}'
                    random_tasks.append(dict(
                        chunk_mtx=chunk_mtx, genes_path=scratch / 'gene_names.txt',
                        features_path=batch_features_path, metadata_path=chunk_meta_path,
                        random_dir=batch_gmt, out_prefix=batch_prefix, include_random=True))
                    batch_outputs.append((batch_prefix, feature_start, feature_count))
                random_specs.append((n_chunk, batch_outputs))
            _assert_wt_random_chunk_nucleus_order(
                metadata.nucleus_id.iloc[wt_row_ids].astype(str).tolist(), wt_random_nucleus_ids)
            _run_ucell_chunks_parallel(random_tasks, args.workers)
            wt_cursor = 0
            for n_chunk, batch_outputs in random_specs:
                for score_prefix, feature_start, feature_count in batch_outputs:
                    for rank in MAX_RANKS:
                        raw_path = Path(f'{score_prefix}.sharecount_{rank}.bin')
                        share_count = np.fromfile(raw_path, dtype=np.uint8)
                        assert share_count.size == n_chunk * feature_count
                        share_count = share_count.reshape((n_chunk, feature_count), order='F')
                        share_memmaps[rank][wt_cursor:wt_cursor + n_chunk,
                                             feature_start:feature_start + feature_count] = share_count
                wt_cursor += n_chunk
            if wt_cursor != len(wt_row_ids) or not len(wt_row_ids):
                raise ValueError(f'{pop}: WT random chunks do not cover every WT discovery nucleus')
            mouse_median_shares_by_rank = {}
            wt_meta = metadata.iloc[wt_row_ids].reset_index(drop=True)
            wt_mouse_labels = wt_meta.mouse_id.astype(str).to_numpy()
            for rank in MAX_RANKS:
                share_memmaps[rank].flush()
                mouse_median = np.empty((len(mouse_ids), len(records)), dtype=np.float64)
                for mouse_i, mouse in enumerate(mouse_ids):
                    rows = np.flatnonzero(wt_mouse_labels == mouse)
                    if not len(rows):
                        raise ValueError(f'{pop}: no WT discovery nuclei in per-nucleus random share matrix for {mouse}')
                    mouse_median[mouse_i] = np.median(np.asarray(share_memmaps[rank][rows], dtype=np.float64), axis=0) / N_RANDOM_SETS
                mouse_median_shares_by_rank[rank] = mouse_median
                share_memmaps[rank].flush()
            score_matrix.flush()
            qc = _summarize_scores(pop, score_matrix, metadata, records, ulm_values,
                                   mouse_median_shares_by_rank, args.out)
            population_qc[pop] = qc['qc_summary']
            population_analysis[pop] = qc
            ulm_estimate.index.name = 'mouse_id'

        bergmann_meta_path = pd.read_csv(args.out / 'state_P_Bergmann_glia_nuclei.tsv', sep='\t',
                                         dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str})
        bergmann_scores = np.load(args.out / 'state_P_Bergmann_glia.npy', mmap_mode='r', allow_pickle=False)
        neighbour_table = pd.read_csv(scratch / 'staged_counts/BERGMANN_T6B_NEIGHBOUR_SCORE.tsv.gz',
                                     sep='\t', dtype={'nucleus_id': str, 'mouse_id': str})
        qc6_summary = run_qc6_bergmann(
            bergmann_scores, bergmann_meta_path, records,
            population_analysis['Bergmann glia']['enriched'][PRIMARY_MAX_RANK],
            neighbour_table, t6b_marker_table,
            scratch / 'staged_counts/Bergmann_glia', scratch / 'gene_names.txt',
            scratch, args.out, workers=args.workers)
        population_qc['Bergmann glia']['QC6'] = qc6_summary
        if qc6_summary['n_flagged_sets']:
            population_qc['Bergmann glia']['qc3_to_qc6_flags'].append(
                f"QC6: {qc6_summary['n_flagged_sets']} enriched set(s) triggered at least one warning condition (LOO-driver marker overlap OR neighbour-score correlation)")

        # Granule is scored only as a leak check; it never contributes to enrichment or ULM.
        granule_features_path = args.out / 'state_P_Purkinje_cell_features.tsv'
        granule_meta_path = scratch / 'staged_counts/granule_leak_only_nuclei.tsv'
        granule_meta = pd.read_csv(granule_meta_path, sep='\t', dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str, 'sex': str})
        granule_sums: dict[str, np.ndarray] = {}
        granule_counts: dict[str, int] = {}
        granule_age: dict[str, int] = {}
        granule_chunks = sorted((scratch / 'staged_counts/granule_leak_only').glob('chunk_*.mtx'))
        granule_summary_dir = scratch / 'granule_mouse_summary'
        granule_summary_dir.mkdir()
        granule_tasks = []
        granule_specs = []
        for granule_chunk_i, chunk_mtx in enumerate(granule_chunks):
            prefix = chunk_mtx.with_suffix('')
            chunk_meta_path = Path(str(prefix) + '.nuclei.tsv')
            score_prefix = granule_summary_dir / f'granule_{granule_chunk_i:06d}'
            granule_tasks.append(dict(
                chunk_mtx=chunk_mtx, genes_path=scratch / 'gene_names.txt',
                features_path=granule_features_path, metadata_path=chunk_meta_path,
                random_dir=None, out_prefix=score_prefix, include_random=False,
                max_ranks=(PRIMARY_MAX_RANK,), mode='mouse_summary'))
            granule_specs.append(score_prefix)
        _run_ucell_chunks_parallel(granule_tasks, args.workers)
        granule_feature_count = len(population_features[POPULATIONS[0]])
        for score_prefix in granule_specs:
            summary_metadata = pd.read_csv(
                Path(f'{score_prefix}.mouse_summary_metadata.tsv'), sep='\t',
                dtype={'mouse_id': str})
            summary_raw = np.fromfile(Path(f'{score_prefix}.mouse_summary.bin'), dtype='<f8')
            if summary_raw.size != len(summary_metadata) * granule_feature_count:
                raise ValueError('Granule mouse-summary score binary has the wrong dimensions')
            summary_matrix = summary_raw.reshape((len(summary_metadata), granule_feature_count), order='F')
            if not np.isfinite(summary_matrix).all():
                raise ValueError('Granule mouse-summary scores contain non-finite values')
            for summary_i, row in summary_metadata.iterrows():
                mouse = str(row.mouse_id)
                age = int(row.age_weeks)
                if mouse in granule_age and granule_age[mouse] != age:
                    raise ValueError(f'Granule mouse-summary age changed between chunks for {mouse}')
                granule_age[mouse] = age
                granule_sums.setdefault(mouse, np.zeros(granule_feature_count, dtype=np.float64))
                granule_sums[mouse] += summary_matrix[summary_i]
                granule_counts[mouse] = granule_counts.get(mouse, 0) + int(row.n_nuclei)
        granule_records = []
        if granule_sums:
            ftab = population_features[POPULATIONS[0]]
            for mouse in sorted(granule_sums):
                means = granule_sums[mouse] / granule_counts[mouse]
                for j, feature in ftab.iterrows():
                    granule_records.append(dict(population=GRANULE, use='leak check only',
                        mouse_id=mouse, age_weeks=granule_age[mouse],
                        feature_id=str(feature.stable_id), set_name=str(feature.set_name),
                        n_nuclei=granule_counts[mouse], mean_UCell_1500=float(means[j])))
        pd.DataFrame(granule_records).to_csv(args.out / 'GRANULE_LEAK_CHECK.tsv.gz', sep='\t', index=False, compression='gzip')

        # Compare qualifying lists at the two sensitivity ranks with the primary list.
        jaccard_all = []
        for pop in POPULATIONS:
            # The per-population rank Jaccards are emitted by _summarize_scores.
            path = args.out / f'{POP_SLUG[pop]}_ENRICHED_SETS.tsv'
            enriched_table = pd.read_csv(path, sep='\t')
            for rank in MAX_RANKS:
                jaccard_all.append(dict(population=pop, maxRank=rank,
                    n_enriched=int((enriched_table.maxRank == rank).sum())))
        pd.DataFrame(jaccard_all).to_csv(args.out / 'ENRICHED_LIST_COUNTS_BY_RANK.tsv', sep='\t', index=False)

        retained_scratch_files = [path for path in scratch.rglob('*') if path.is_file()]
        retained_scratch_bytes = sum(path.stat().st_size for path in retained_scratch_files)
        expected_fixed_ucell_scratch_bytes = sum(
            int(scoring_plan[pop]['expected_retained_UCell_binary_bytes']) for pop in POPULATIONS)
        staged_count_bytes = tree_file_size_bytes(scratch / 'staged_counts')
        random_control_gmt_bytes = sum(
            tree_file_size_bytes(scratch / POP_SLUG[pop] / 'random_gmt') for pop in POPULATIONS)
        granule_mouse_summary_bytes = tree_file_size_bytes(scratch / 'granule_mouse_summary')
        expected_core_scratch_bytes = (staged_count_bytes + random_control_gmt_bytes
                                       + expected_fixed_ucell_scratch_bytes
                                       + granule_mouse_summary_bytes)
        expected_scratch_size = dict(
            staged_count_and_metadata_bytes=staged_count_bytes,
            random_control_GMT_bytes=random_control_gmt_bytes,
            fixed_UCell_binary_bytes=expected_fixed_ucell_scratch_bytes,
            Granule_per_mouse_summary_bytes=granule_mouse_summary_bytes,
            expected_core_bytes=expected_core_scratch_bytes,
            headroom_fraction=0.25,
            expected_with_headroom_bytes=(expected_core_scratch_bytes * 5 + 3) // 4,
            headroom_for=['QC6 LOO intermediates', 'logs', 'small metadata and text outputs',
                          'gene-name index and bookkeeping files'],
            all_scratch_files_retained=True)
        analysis_disposition = ('STOP' if any(q['hard_gates']['hard_gate_status'] == 'STOP'
                                               for q in population_qc.values()) else 'PASS')
        record = dict(mode='real', scope='gene_set_space_P', completed=True,
            analysis_disposition=analysis_disposition,
            input_sha256=pins,
            nucleus_metadata_source=dict(selection='included discovery nuclei at 5/18/24/30 plus the six SCOPING-locked 12-week mice, age-only re-inclusion',
                nucleus_manifest_sha256=pins[str(args.manifest)],
                animal_manifest_sha256=pins[str(args.animal_manifest)],
                scoping_record_sha256=pins[str(args.scoping)],
                metadata_counts='log10_total_umi and n_detected from raw full layers/counts rows after sum_duplicates and before symbol aggregation',
                per_population_12_week_counts=TWELVE_POPULATION_COUNTS),
            environment=env_record,
            R_runtime_probe=dict(invocation='Rscript --vanilla',
                libPaths=env_record['R_libPaths'], packageVersion_UCell=env_record['R_package_versions']['UCell']),
            script_sha256=script_hashes(), approval_gate=gate, worker_count=args.workers,
            random_chunk_rows=args.random_chunk_rows,
            execution_backend='local subprocess',
            python_launch=dict(invocation='python -I -B', isolated=sys.flags.isolated == 1,
                               dont_write_bytecode=bool(sys.dont_write_bytecode)),
            T6B_marker_definition=dict(script_path=str(T6B_SCRIPT), script_sha256=T6B_SCRIPT_SHA256,
                run_folder_read=False, marker_table='T6B_WT_TOP_MARKERS.tsv',
                marker_score_file='T6B_MARKER_DEFINITION.json'),
            duplicate_symbol_aggregation={k: source_details[k] for k in (
                'n_original_gene_columns', 'n_unique_gene_symbols',
                'n_duplicate_symbol_columns_aggregated', 'n_distinct_duplicated_symbols',
                'duplicate_symbol_names')},
            scoring_plan=scoring_plan,
            scratch_retention=dict(path=str(scratch), files_retained=True,
                retained_file_count=len(retained_scratch_files), retained_bytes=retained_scratch_bytes,
                expected_fixed_UCell_binary_bytes=expected_fixed_ucell_scratch_bytes,
                expected_fixed_bytes_excludes=['staged MatrixMarket counts', 'random GMT controls',
                                               'QC6 LOO effect/base-score files', 'text logs and metadata']),
            expected_scratch_size=expected_scratch_size,
            QC_hard_gates={pop: qc['hard_gates'] for pop, qc in population_qc.items()},
            ULM_missing_by_population=ulm_missing_by_population,
            QC6_Bergmann=qc6_summary)
        _write_json(args.out / 'RECORD.json', record)
        if analysis_disposition == 'STOP':
            _write_json(args.out / 'STOP.json', dict(
                status='STOP', reason='QC1 or QC2 hard gate failed; outputs are preserved for reporting, not interpretation.',
                QC_hard_gates={pop: qc['hard_gates'] for pop, qc in population_qc.items()}))

        # RESULT carries the approved settings, analysis disposition and legacy-UCell disclosure.
        result = dict(mode='real', scope='gene_set_space_P', completed=True,
            analysis_disposition=analysis_disposition,
            input_sha256=pins,
            nucleus_metadata_source=record['nucleus_metadata_source'],
            worker_count=args.workers, random_chunk_rows=args.random_chunk_rows,
            execution_backend='local subprocess',
            scoring_plan=scoring_plan,
            scratch_retention=record['scratch_retention'],
            expected_scratch_size=expected_scratch_size,
            gene_set_source='native mouse MSigDB 2025.1.Mm via msigdbr 25.1.1; no human ortholog conversion',
            gmt_sha256=GMT_SHA256, official_gmt_comparison=comparison,
            gene_symbol_handling=record['duplicate_symbol_aggregation'],
            populations=list(POPULATIONS), granule_role='leak check only; excluded from enrichment and ULM',
            candidate_rule=dict(min_matched_genes=MIN_MATCHED_GENES, max_matched_genes=MAX_MATCHED_GENES,
                                missing_genes='unmatched native symbols are not imputed as expressed; the matched-gene set is scored and absent assay genes have zero counts'),
            UCell=dict(version=UCELL_VERSION, primary_maxRank=PRIMARY_MAX_RANK, sweep_maxRanks=[1000, 3000],
                      chunk_size=CHUNK_ROWS, WT_random_chunk_rows=args.random_chunk_rows, seed=SEED, smoothing='none',
                      disclosure='UCell 2.6.2 predates the 2026 “UCell v2”, which corrected the theoretical maximum of U. This changes only the scale within a set, not the ordering of nuclei, and the random-set null uses sets of the same size.'),
            python_launch=dict(invocation='python -I -B', isolated=sys.flags.isolated == 1,
                               dont_write_bytecode=bool(sys.dont_write_bytecode)),
            random_null=dict(n_random_sets=N_RANDOM_SETS, expression_bins=EXPRESSION_BINS,
                             seed=SEED, size_matched=True, expression_bin_matched=True,
                             top_n_cap=None,
                             statistic='For each nucleus and set, fraction of 200 matched random-set UCell scores strictly beaten by the observed score; the mouse statistic is the median across its WT discovery nuclei.',
                             criterion='mouse median nucleus-share >=0.95 in at least 11/12 WT discovery mice, AND positive ULM in those same mice'),
            ULM=dict(method='decoupler run_ulm', version=DECOUPLER_VERSION,
                     min_n=MIN_MATCHED_GENES,
                     missing_source_policy='A source absent after decoupler min_n filtering is recorded as ULM-missing and assigned zero ULM values before enrichment; it cannot be ULM-positive.',
                     input='per-mouse mean of nucleus-level log1p(CP10k), WT discovery ages 5/18/24/30 only'),
            ULM_missing_by_population=ulm_missing_by_population,
            QC_reporting_flags=dict(
                QC1='STOP if median across nuclei of set count with UCell score >=0.99 is >0',
                QC2='STOP if >1% of nuclei have within-nucleus UCell score IQR exactly 0',
                QC3='flag any set with |Spearman rho|>0.4 against log10 total UMI or detected genes',
                QC4=f'flag any fixed-seed 50% WT downsample mouse-mean Pearson r <{QC4_DOWNSAMPLE_R_FLAG} or undefined',
                QC5=f'flag any pairwise enriched-list maxRank Jaccard <{QC5_MAXRANK_JACCARD_FLAG}',
                QC6='flag a Bergmann enriched set if either a top-decile leave-one-gene-out driver overlaps the T6b top Granule/Purkinje marker list OR WT Bergmann set-score Spearman rho is >=0.4 against the T6b neighbour score; record which condition fired',
                QC3_to_QC6_are_reporting_flags_not_hard_gates=True,
                exact_ceiling_and_zero_shares='reported separately for every matched set'),
            QC6_Bergmann=dict(**qc6_summary,
                T6b_script_sha256=T6B_SCRIPT_SHA256,
                Granule_marker_ratio_threshold=T6B_GRANULE_QC6_RATIO,
                Granule_marker_count=granule_marker_count,
                Granule_marker_count_at_ratio_4=granule_marker_count_at_ratio4,
                Purkinje_marker_ratio_threshold=T6B_PRIMARY_RATIO,
                Purkinje_marker_count=purkinje_marker_count,
                Bergmann_self_marker_ratio_threshold=T6B_PRIMARY_RATIO,
                Bergmann_self_marker_count=bergmann_marker_count,
                Granule_marker_count_at_ratio_2=49,
                Granule_threshold_reason='QC6-only predeclared relaxation: the pinned T6b tables show 0 Granule markers at ratio 4, 6 at ratio 3, and 49 at ratio 2; ratio 3 was selected as the first relaxation yielding the predeclared six markers and is used only for this warning check.'),
            QC=population_qc,
            environment=env_record, script_sha256=script_hashes(), approval_gate=gate,
            RECORD='RECORD.json',
            result_files=sorted(p.name for p in args.out.iterdir() if p.is_file()
                                and p.name not in {'RESULT.json', 'RECORD.json'}),
            elapsed_seconds=time.time() - run_started)
        result['official_gmt_comparison']['pinned_input_remains'] = '712fd1c0f9dcec85b45da2aeb75a355d8f3c22542f0b1a77cd4520d89d2b222d'
        _write_json(args.out / 'RESULT.json', result)
        (args.out / 'CONDA_EXPLICIT.txt').write_bytes(subprocess.run(
            [str(CONDA), 'list', '--explicit', '-p', str(ENV)], check=True, capture_output=True).stdout)
        write_p_interface_checksum_manifest(args.out)
        return result
    except BaseException:
        _write_json(args.out / 'FAILED.json', dict(mode='real', scope='gene_set_space_P',
            script_sha256=script_hashes(), traceback=traceback.format_exc(),
            single_shot='failed output is preserved; do not retry or overwrite'))
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--real', action='store_true', required=True)
    parser.add_argument('--go-file', type=Path, required=True)
    parser.add_argument('--h5ad', type=Path, default=H5AD)
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--animal-manifest', type=Path, default=ANIMAL_MANIFEST)
    parser.add_argument('--scoping', type=Path, default=SCOPING)
    parser.add_argument('--gmt', type=Path, default=GMT)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--workers', type=int, default=DEFAULT_UCELL_PARALLEL_WORKERS,
                        help='maximum concurrent Rscript tasks (default: 4)')
    parser.add_argument('--random-chunk-rows', type=int, default=DEFAULT_RANDOM_CHUNK_ROWS,
                        help='WT-discovery-only random-scoring chunk rows (default: 256)')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be a positive integer')
    if not 1 <= args.random_chunk_rows <= CHUNK_ROWS:
        parser.error(f'--random-chunk-rows must be between 1 and {CHUNK_ROWS}')
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('refusing: critical checks require Python assertions')
    if sys.flags.isolated != 1 or not sys.dont_write_bytecode:
        raise SystemExit('refusing real-data entry point unless launched with python -I -B')
    try:
        result = run_real(args)
        print(json.dumps(dict(completed=True, output=str(args.out), files=len(result['result_files'])), sort_keys=True))
    except BaseException:
        raise


if __name__ == '__main__':
    main()
