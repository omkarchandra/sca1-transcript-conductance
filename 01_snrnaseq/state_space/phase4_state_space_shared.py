"""Phase 4 shared P/G path, mouse readout, and WT-only DPT module.

Draft only. Real inputs are opened only after a shared_readout-scoped gate binds
this script hash, the exact P/G interface checksum files, review byte review,
and <AUTHOR>'s run-specific literal written ``yes``. Input matrices use
state_<space>_<celltype>.npy with a matching nuclei.tsv and features.tsv.
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
import importlib.metadata
import itertools
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
POPS = ('Purkinje cell', 'Bergmann glia')
POP_SLUG = {p: p.replace(' ', '_') for p in POPS}
DISCOVERY_AGES = (5, 18, 24, 30)
READOUT_AGES = (5, 12, 18, 24, 30)
N_PCS = 20
N_NEIGHBORS = 15
N_DCS = 10
DPT_SEED = 20261008
PCA_SEED = 20261008
PCA_DIAGNOSTIC_SEED = 20261008
PCA_DIAGNOSTIC_NUCLEI = 2000
CHUNK = 256
SPEARMAN_FLOOR = 0.6
ENV = Path(_release_paths.expand('${CONDA_ENV}'))
CONDA = Path(_release_paths.expand('${CONDA_BIN}'))
CONDA_EXPLICIT_SHA256 = '83a1d457bc19a6e80ad7613bc16ba83a01b8f92a73ab992e8185161950d703a0'
CONDA_EXPLICIT_BYTES = 43839
ENV_VERSIONS = {'python': '3.9.23', 'scanpy': '1.10.3', 'numpy': '1.26.4', 'scipy': '1.13.1'}
BUILD_DOCS = ROOT / _release_paths.expand('${DATA_ROOT}/phase4_state_space_build')
REAL_OUT = Path(_release_paths.expand('${WORK_ROOT}/phase4_state_space_shared'))


def sha(path: Union[Path, str]) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(16 << 20), b''):
            h.update(block)
    return h.hexdigest()


def json_safe(value):
    """Convert NumPy/Pandas scalars and nonfinite numbers to strict JSON values."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def check_go(path: Union[Path, str], expected_scope: str,
             expected_scripts: dict[str, str], p_dir: Optional[Path] = None,
             g_dir: Optional[Path] = None) -> dict:
    """Validate one scoped approval; shared runs bind both interface checksum files."""
    gate = json.loads(Path(path).read_text(encoding='utf-8'))
    assert expected_scope in {'G_builder', 'shared_readout'}
    assert gate.get('scope') == expected_scope, ('gate scope mismatch', gate.get('scope'), expected_scope)
    assert gate.get('byte_clearance') == 'yes', 'byte review is absent'
    assert gate.get('author_written_yes') == 'yes', '<AUTHOR> literal written yes is absent'
    assert gate.get('script_sha256') == expected_scripts, ('gate is bound to different script bytes', gate.get('script_sha256'), expected_scripts)
    assert gate.get('relayed_by') and gate.get('relayed_at_cdt'), 'approval relay provenance is missing'
    if expected_scope == 'shared_readout':
        assert p_dir is not None and g_dir is not None, 'shared-readout gate check needs --p-dir and --g-dir'
        p_sums = Path(p_dir) / 'SHA256SUMS.txt'
        g_sums = Path(g_dir) / 'SHA256SUMS.txt'
        assert p_sums.is_file() and g_sums.is_file(), 'an input interface SHA256SUMS.txt is missing'
        p_digest, g_digest = sha(p_sums), sha(g_sums)
        assert gate.get('p_interface_sha256sums') == p_digest, ('P interface checksum file differs from gate', p_digest)
        assert gate.get('g_interface_sha256sums') == g_digest, ('G interface checksum file differs from gate', g_digest)
        gate['verified_interface_sha256sums'] = {
            'p_interface_sha256sums': p_digest,
            'g_interface_sha256sums': g_digest,
        }
    return gate


def validate_env_provenance(capture: bool = False) -> tuple[dict, Optional[bytes]]:
    """Check the unchanged cluster environment's explicit export and core versions."""
    command = [str(CONDA), 'list', '--explicit', '-p', str(ENV)]
    exported = subprocess.run(command, check=True, capture_output=True).stdout
    got_sha = hashlib.sha256(exported).hexdigest()
    assert len(exported) == CONDA_EXPLICIT_BYTES and got_sha == CONDA_EXPLICIT_SHA256, ('single_cell explicit export changed', got_sha)
    versions = dict(python=sys.version.split()[0], numpy=np.__version__)
    assert versions['python'] == ENV_VERSIONS['python'] and versions['numpy'] == ENV_VERSIONS['numpy']
    assert pd.__version__ == '2.3.3'
    # Scanpy is imported lazily after the command is launched with the env's lib path.
    import scanpy as sc
    import scipy
    versions.update(scanpy=importlib.metadata.version('scanpy'), scipy=scipy.__version__)
    assert versions == ENV_VERSIONS, ('core package versions changed', versions, ENV_VERSIONS)
    record = dict(environment_path=str(ENV), export_command='conda list --explicit -p ' + str(ENV),
                  explicit_file='CONDA_EXPLICIT.txt', explicit_sha256=got_sha, explicit_bytes=len(exported),
                  versions=versions, no_environment_changes=True)
    return record, exported if capture else None


def _mouse_weights(meta: pd.DataFrame, rows: np.ndarray) -> np.ndarray:
    mice = meta.iloc[rows].mouse_id.astype(str).to_numpy()
    unique, inverse, counts = np.unique(mice, return_inverse=True, return_counts=True)
    assert len(unique) >= 2 and np.all(counts > 0)
    return 1.0 / (len(unique) * counts[inverse])


def _z_chunk(X: np.ndarray, rows: np.ndarray, model: dict) -> np.ndarray:
    raw = np.asarray(X[rows], dtype=np.float64)
    umi = model['umi'][rows]
    residual = raw - model['intercept'][None, :] - umi[:, None] * model['slope'][None, :]
    return np.clip((residual - model['residual_mean'][None, :]) / model['residual_sd'][None, :], -10.0, 10.0)


def fit_wt_model(X: np.ndarray, meta: pd.DataFrame, train_mask: np.ndarray, n_pcs: int = N_PCS,
                 seed: int = PCA_SEED) -> dict:
    """Mouse-balanced UMI regression, z-score/clip, then deterministic weighted randomized PCA."""
    train_mask = np.asarray(train_mask, bool)
    rows = np.flatnonzero(train_mask)
    assert len(rows) and X.ndim == 2 and X.shape[0] == len(meta)
    assert meta.iloc[rows].genotype.eq('WT').all(), 'KI nuclei entered a WT fit'
    assert meta.iloc[rows].age_weeks.isin(DISCOVERY_AGES).all(), '12 W or other ages entered the fit'
    weights = _mouse_weights(meta, rows)
    p = X.shape[1]
    umi_all = meta.log10_total_umi.to_numpy(dtype=np.float64)
    umi = umi_all[rows]
    # The requested UMI regression and feature z-scores are nucleus-level;
    # equal-mouse weights are applied specifically to the PCA covariance.
    xbar = float(np.mean(umi))
    xvar = float(np.mean((umi - xbar) ** 2))
    assert np.isfinite(xvar) and xvar > 0, 'log10 total UMI has no training variation'
    sum_y = np.zeros(p, dtype=np.float64)
    sum_xy = np.zeros(p, dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        rr = rows[first:first + CHUNK]
        values = np.asarray(X[rr], dtype=np.float64)
        sum_y += values.sum(axis=0)
        sum_xy += values.T @ umi[first:first + len(rr)]
    ybar = sum_y / len(rows)
    sum_xy /= len(rows)
    slope = (sum_xy - xbar * ybar) / xvar
    intercept = ybar - slope * xbar
    sum_residual = np.zeros(p, dtype=np.float64)
    sum_residual2 = np.zeros(p, dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        rr = rows[first:first + CHUNK]
        values = np.asarray(X[rr], dtype=np.float64)
        residual = values - intercept[None, :] - umi[first:first + len(rr), None] * slope[None, :]
        sum_residual += residual.sum(axis=0)
        sum_residual2 += (residual * residual).sum(axis=0)
    residual_mean = sum_residual / len(rows)
    residual_sd = np.sqrt(np.maximum(sum_residual2 / len(rows) - residual_mean ** 2, 0.0))
    residual_sd[residual_sd == 0] = 1.0
    model = dict(umi=umi_all, slope=slope, intercept=intercept,
                 residual_mean=residual_mean, residual_sd=residual_sd)
    sum_z = np.zeros(p, dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        rr = rows[first:first + CHUNK]
        w = weights[first:first + len(rr)]
        sum_z += _z_chunk(X, rr, model).T @ w
    z_mean = sum_z
    model['z_mean'] = z_mean
    total_variance = 0.0
    for first in range(0, len(rows), CHUNK):
        rr = rows[first:first + CHUNK]
        z = _z_chunk(X, rr, model) - z_mean
        w = weights[first:first + len(rr)]
        total_variance += float(np.sum(w * np.einsum('ij,ij->i', z, z)))
    assert total_variance > 0 and np.isfinite(total_variance)
    k = min(int(n_pcs), p)
    assert k > 0 and len(rows) > k, 'insufficient nuclei or features for requested PCA rank'
    ell = min(len(rows), p, k + min(10, max(2, k // 2)))
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((p, ell))
    sqrt_weights = np.sqrt(weights)
    Y = np.empty((len(rows), ell), dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        last = min(first + CHUNK, len(rows)); rr = rows[first:last]
        z = _z_chunk(X, rr, model) - z_mean
        Y[first:last] = (z * sqrt_weights[first:last, None]) @ omega
    Q, _ = np.linalg.qr(Y, mode='reduced')
    # Two deterministic subspace iterations improve recovery of lower retained PCs.
    for _ in range(2):
        right = np.zeros((p, Q.shape[1]), dtype=np.float64)
        for first in range(0, len(rows), CHUNK):
            last = min(first + CHUNK, len(rows)); rr = rows[first:last]
            z = _z_chunk(X, rr, model) - z_mean
            right += z.T @ (sqrt_weights[first:last, None] * Q[first:last])
        next_y = np.empty_like(Q)
        for first in range(0, len(rows), CHUNK):
            last = min(first + CHUNK, len(rows)); rr = rows[first:last]
            z = _z_chunk(X, rr, model) - z_mean
            next_y[first:last] = (z @ right) * sqrt_weights[first:last, None]
        Q, _ = np.linalg.qr(next_y, mode='reduced')
    B = np.zeros((Q.shape[1], p), dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        last = min(first + CHUNK, len(rows)); rr = rows[first:last]
        z = _z_chunk(X, rr, model) - z_mean
        B += Q[first:last].T @ (z * sqrt_weights[first:last, None])
    _u, singular, vt = np.linalg.svd(B, full_matrices=False)
    rotation = vt[:k].T.copy()
    eigenvalues = singular[:k] ** 2
    variance_fractions = eigenvalues / total_variance
    for j in range(k):
        anchor = int(np.argmax(np.abs(rotation[:, j])))
        if rotation[anchor, j] < 0:
            rotation[:, j] *= -1
    model.update(rotation=rotation, eigenvalues=eigenvalues,
                 variance_fractions=variance_fractions,
                 cumulative_variance_fractions=np.cumsum(variance_fractions),
                 total_weighted_variance=total_variance,
                 n_pcs=k, n_features=p, n_fit_nuclei=len(rows),
                 n_fit_mice=int(meta.iloc[rows].mouse_id.nunique()),
                 pca_seed=int(seed), pca_power_iterations=2,
                 pca_oversamples=int(ell - k),
                 pca_method='mouse-weighted randomized SVD; 2 power iterations; up to 10 oversamples')
    return model


def pca_diagnostics(X: np.ndarray, meta: pd.DataFrame, train_mask: np.ndarray,
                    full_model: dict, diagnostic_n: int = PCA_DIAGNOSTIC_NUCLEI) -> tuple[pd.DataFrame, dict]:
    """Report full-fit variance and compare randomized/exact loadings on a fixed WT subsample."""
    rows = np.flatnonzero(np.asarray(train_mask, dtype=bool))
    requested_n = int(diagnostic_n)
    assert requested_n == PCA_DIAGNOSTIC_NUCLEI or requested_n > full_model['n_pcs'], (
        'diagnostic subsample must be the declared production size or exceed the retained rank')
    sample_n = min(requested_n, len(rows))
    used_all_wt_discovery_nuclei = bool(sample_n == len(rows))
    assert sample_n > full_model['n_pcs'], 'insufficient WT-discovery nuclei for exact PCA loading check'
    assert sample_n > full_model['n_pcs']
    rng = np.random.default_rng(PCA_DIAGNOSTIC_SEED)
    sample_rows = np.sort(rng.choice(rows, size=sample_n, replace=False))
    sample_mask = np.zeros(len(meta), dtype=bool)
    sample_mask[sample_rows] = True
    sample_model = fit_wt_model(X, meta, sample_mask, n_pcs=N_PCS, seed=PCA_SEED)
    weights = _mouse_weights(meta, sample_rows)
    weighted_z = np.empty((sample_n, X.shape[1]), dtype=np.float64)
    for first in range(0, sample_n, CHUNK):
        last = min(first + CHUNK, sample_n)
        rr = sample_rows[first:last]
        z = _z_chunk(X, rr, sample_model) - sample_model['z_mean']
        weighted_z[first:last] = z * np.sqrt(weights[first:last, None])
    _u, _singular, exact_vt = np.linalg.svd(weighted_z, full_matrices=False)
    n_compare = min(5, sample_model['n_pcs'], exact_vt.shape[0])
    cosines = np.clip(
        np.abs(np.sum(sample_model['rotation'][:, :n_compare] * exact_vt[:n_compare].T, axis=0)),
        0.0, 1.0)
    rows_out = []
    for pc in range(full_model['n_pcs']):
        rows_out.append(dict(
            diagnostic='full_fit_variance', pc=pc + 1,
            eigenvalue=float(full_model['eigenvalues'][pc]),
            variance_fraction=float(full_model['variance_fractions'][pc]),
            cumulative_variance_fraction=float(full_model['cumulative_variance_fractions'][pc]),
            fixed_seed_subsample_n=sample_n, subsample_seed=PCA_DIAGNOSTIC_SEED,
            used_all_wt_discovery_nuclei=used_all_wt_discovery_nuclei,
            randomized_vs_exact_abs_cosine=float(cosines[pc]) if pc < n_compare else np.nan))
    nucleus_ids = meta.iloc[sample_rows].nucleus_id.astype(str).tolist()
    result = dict(
        total_weighted_variance=float(full_model['total_weighted_variance']),
        n_pcs=int(full_model['n_pcs']),
        variance_captured_fraction=float(np.sum(full_model['variance_fractions'])),
        fixed_seed_subsample_n=int(sample_n),
        fixed_seed_subsample_requested_n=requested_n,
        used_all_wt_discovery_nuclei=used_all_wt_discovery_nuclei,
        fixed_seed_subsample_seed=int(PCA_DIAGNOSTIC_SEED),
        fixed_seed_subsample_nucleus_ids_sha256=hashlib.sha256('\n'.join(nucleus_ids).encode()).hexdigest(),
        exact_svd_scope='sampled WT nuclei only',
        exact_loading_cosines={f'PC{i + 1}': float(cosines[i]) for i in range(n_compare)},
        n_exact_loading_cosines=int(n_compare))
    return pd.DataFrame(rows_out), result


def transform_pcs(X: np.ndarray, rows: np.ndarray, model: dict) -> np.ndarray:
    rows = np.asarray(rows, dtype=int)
    out = np.empty((len(rows), model['n_pcs']), dtype=np.float64)
    for first in range(0, len(rows), CHUNK):
        rr = rows[first:first + CHUNK]
        z = _z_chunk(X, rr, model) - model['z_mean']
        out[first:first + len(rr)] = z @ model['rotation']
    return out


def _mouse_medians(points: np.ndarray, meta: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    frame = pd.DataFrame(points, columns=[f'PC{i + 1}' for i in range(points.shape[1])])
    frame['mouse_id'] = meta.mouse_id.astype(str).to_numpy()
    frame['age_weeks'] = meta.age_weeks.to_numpy(dtype=int)
    frame['genotype'] = meta.genotype.to_numpy()
    frame['nucleus_id'] = meta.nucleus_id.astype(str).to_numpy()
    med = frame.groupby(['mouse_id', 'age_weeks', 'genotype'], sort=True)[[f'PC{i + 1}' for i in range(points.shape[1])]].median().reset_index()
    return med[[f'PC{i + 1}' for i in range(points.shape[1])]].to_numpy(), med


def age_centers(points: np.ndarray, meta: pd.DataFrame) -> np.ndarray:
    """Each age centre is the mean of the WT per-mouse medians."""
    _values, med = _mouse_medians(points, meta.reset_index(drop=True))
    centers = []
    for age in DISCOVERY_AGES:
        g = med[(med.age_weeks == age) & (med.genotype == 'WT')]
        assert len(g) >= 2, f'age {age}: fewer than two WT fitting mice in this fit'
        cols = [f'PC{i + 1}' for i in range(points.shape[1])]
        centers.append(g[cols].to_numpy(dtype=np.float64).mean(axis=0))
    return np.asarray(centers)


def project_to_path(points: np.ndarray, centers: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nearest polyline projection; first/last segments extend as rays in t."""
    points = np.asarray(points, dtype=np.float64)
    centers = np.asarray(centers, dtype=np.float64)
    assert centers.shape == (4, points.shape[1])
    best_sq = np.full(len(points), np.inf)
    best_t = np.full(len(points), np.nan)
    best_foot = np.empty_like(points)
    for j, (age0, age1) in enumerate(zip(DISCOVERY_AGES[:-1], DISCOVERY_AGES[1:])):
        a = centers[j]; v = centers[j + 1] - a
        vv = float(np.dot(v, v))
        assert vv > 0, f'zero-length age-path segment {age0}-{age1}'
        alpha = ((points - a) @ v) / vv
        lo = -np.inf if j == 0 else 0.0
        hi = np.inf if j == len(DISCOVERY_AGES) - 2 else 1.0
        alpha = np.clip(alpha, lo, hi)
        foot = a[None, :] + alpha[:, None] * v[None, :]
        sq = np.einsum('ij,ij->i', points - foot, points - foot)
        take = sq < best_sq - 1e-14
        best_sq[take] = sq[take]
        best_t[take] = age0 + alpha[take] * (age1 - age0)
        best_foot[take] = foot[take]
    return best_t, np.sqrt(best_sq), best_foot


def _scanpy_dpt(graph_points: np.ndarray, graph_meta: pd.DataFrame, seed: int = DPT_SEED) -> tuple[np.ndarray, dict]:
    """Build the WT-only balanced kNN graph and return graph-node Scanpy DPT."""
    assert len(graph_points) > max(N_DCS + 1, N_NEIGHBORS)
    assert graph_meta.genotype.eq('WT').all(), 'non-WT nucleus entered the DPT graph'
    assert graph_meta.age_weeks.isin(DISCOVERY_AGES).all(), 'non-discovery nucleus entered the DPT graph'
    import scanpy as sc
    import anndata as ad
    from scipy.sparse.csgraph import connected_components

    graph_meta = graph_meta.reset_index(drop=True)
    counts = graph_meta.groupby('mouse_id').size()
    per_mouse = int(counts.min())
    rng = np.random.default_rng(seed)
    chosen: list[int] = []
    ids = graph_meta.nucleus_id.astype(str).to_numpy()
    # Stable ID order makes the fixed-seed sample independent of input row order.
    for mouse in sorted(counts.index.astype(str)):
        candidates = np.flatnonzero(graph_meta.mouse_id.astype(str).to_numpy() == mouse)
        candidates = candidates[np.argsort(ids[candidates], kind='stable')]
        if len(candidates) == per_mouse:
            take = candidates
        else:
            take = rng.choice(candidates, size=per_mouse, replace=False)
        chosen.extend(map(int, take))
    graph_rows = np.asarray(sorted(chosen, key=lambda i: (str(graph_meta.iloc[i].mouse_id), ids[i])), dtype=int)
    sub_meta = graph_meta.iloc[graph_rows].reset_index(drop=True)
    sub_points = np.asarray(graph_points, dtype=np.float64)[graph_rows]
    training_centers = age_centers(np.asarray(graph_points, dtype=np.float64), graph_meta)
    age5 = np.flatnonzero(sub_meta.age_weeks.to_numpy(dtype=int) == 5)
    assert len(age5), 'no WT 5 W DPT graph nuclei'
    center5 = training_centers[0]
    root = min(age5, key=lambda i: (float(np.sum((sub_points[i] - center5) ** 2)), str(sub_meta.iloc[i].nucleus_id)))
    adata = ad.AnnData(X=sub_points.astype(np.float64, copy=False))
    adata.obs_names = sub_meta.nucleus_id.astype(str).to_numpy()
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=0, use_rep='X', method='umap', random_state=int(seed))
    sc.tl.diffmap(adata, n_comps=N_DCS, neighbors_key='neighbors')
    adata.uns['iroot'] = int(root)
    sc.tl.dpt(adata, n_dcs=N_DCS, n_branchings=0, neighbors_key='neighbors')
    dpt = adata.obs['dpt_pseudotime'].to_numpy(dtype=np.float64)
    assert not np.isnan(dpt).any(), 'Scanpy DPT returned NaN; disconnected values must be inf'
    conn = adata.obsp['connectivities']
    n_components, labels = connected_components(conn, directed=False, return_labels=True)
    info = dict(scanpy_version=importlib.metadata.version('scanpy'), n_graph_nodes=len(sub_meta), n_per_mouse=per_mouse,
                n_mice=int(sub_meta.mouse_id.nunique()), k=N_NEIGHBORS, n_dcs=N_DCS,
                random_seed=int(seed), root_nucleus_id=str(sub_meta.iloc[root].nucleus_id),
                root_age_weeks=5, n_connected_components=int(n_components),
                n_graph_infinite=int(np.isinf(dpt).sum()))
    return dpt, dict(graph_rows=graph_rows, graph_points=sub_points, graph_meta=sub_meta,
                     graph_dpt=dpt, info=info, component_labels=labels)


def _transfer_dpt(query_points: np.ndarray, graph_points: np.ndarray, graph_dpt: np.ndarray,
                  n_neighbors: int = N_NEIGHBORS) -> np.ndarray:
    assert len(graph_points) >= n_neighbors
    model = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto', metric='euclidean')
    model.fit(graph_points)
    output = np.empty(len(query_points), dtype=np.float64)
    for first in range(0, len(query_points), 8192):
        last = min(first + 8192, len(query_points))
        idx = model.kneighbors(query_points[first:last], return_distance=False)
        output[first:last] = np.mean(graph_dpt[idx], axis=1)
    return output


def dpt_from_wt_graph(train_points: np.ndarray, train_meta: pd.DataFrame,
                      query_points: np.ndarray, query_meta: pd.DataFrame,
                      seed: int = DPT_SEED) -> tuple[np.ndarray, dict]:
    """Fit a balanced WT-discovery Scanpy graph and transfer DPT to query rows."""
    train_meta = train_meta.reset_index(drop=True)
    assert train_meta.genotype.eq('WT').all() and train_meta.age_weeks.isin(DISCOVERY_AGES).all()
    assert len(train_points) == len(train_meta) and len(query_points) == len(query_meta)
    _sub_dpt, graph = _scanpy_dpt(np.asarray(train_points, dtype=np.float64), train_meta, seed=seed)
    graph_ids = graph['graph_meta'].nucleus_id.astype(str).to_numpy()
    lookup = {name: i for i, name in enumerate(graph_ids)}
    query_ids = query_meta.nucleus_id.astype(str).to_numpy()
    dpt = np.full(len(query_ids), np.nan, dtype=np.float64)
    direct = np.array([name in lookup for name in query_ids], dtype=bool)
    if direct.any():
        dpt[direct] = np.asarray([graph['graph_dpt'][lookup[name]] for name in query_ids[direct]])
    if (~direct).any():
        dpt[~direct] = _transfer_dpt(np.asarray(query_points, dtype=np.float64)[~direct], graph['graph_points'], graph['graph_dpt'])
    assert not np.isnan(dpt).any()
    info = dict(graph['info'], n_query=len(query_meta), n_query_graph_nodes=int(direct.sum()),
                n_query_transferred=int((~direct).sum()), n_query_infinite=int(np.isinf(dpt).sum()))
    graph_nodes = graph['graph_meta'][['nucleus_id', 'mouse_id', 'age_weeks', 'genotype']].copy()
    graph_nodes['graph_DPT'] = graph['graph_dpt']
    graph_nodes['connected_component'] = graph['component_labels']
    graph_nodes['is_root'] = graph_nodes.nucleus_id.astype(str).eq(info['root_nucleus_id']).to_numpy()
    info['graph_nodes'] = graph_nodes.to_dict(orient='records')
    return dpt, info


def balanced_dpt(points: np.ndarray, meta: pd.DataFrame, fit_mask: np.ndarray, query_points: np.ndarray,
                 query_meta: pd.DataFrame, seed: int = DPT_SEED) -> tuple[np.ndarray, dict]:
    """Array-aligned wrapper for callers that keep all rows in one matrix."""
    fit_rows = np.flatnonzero(np.asarray(fit_mask, bool))
    fit_meta = meta.iloc[fit_rows].reset_index(drop=True)
    return dpt_from_wt_graph(np.asarray(points)[fit_rows], fit_meta,
                             np.asarray(query_points), query_meta, seed=seed)


def mean_knn_distance(query: np.ndarray, reference: np.ndarray, k: int = N_NEIGHBORS) -> np.ndarray:
    assert len(reference) >= k, f'need {k} reference nuclei, found {len(reference)}'
    model = NearestNeighbors(n_neighbors=k, metric='euclidean', algorithm='auto').fit(reference)
    out = np.empty(len(query), dtype=np.float64)
    for first in range(0, len(query), 8192):
        last = min(first + 8192, len(query))
        distances, _ = model.kneighbors(query[first:last], return_distance=True)
        out[first:last] = distances.mean(axis=1)
    return out


def leave_one_wt_mouse_reference_rows(meta: pd.DataFrame, age: int) -> list[tuple[str, np.ndarray]]:
    age = int(age)
    age_wt = meta.age_weeks.astype(int).eq(age).to_numpy() & meta.genotype.eq('WT').to_numpy()
    mice = sorted(meta.loc[age_wt, 'mouse_id'].astype(str).unique())
    assert len(mice) == 3, (age, mice, 'expected three same-age WT references')
    mouse_ids = meta.mouse_id.astype(str).to_numpy()
    return [(heldout, np.flatnonzero(age_wt & (mouse_ids != heldout))) for heldout in mice]


def mean_knn_distance_equal_loo_references(query_points: np.ndarray, full_points: np.ndarray,
                                           meta: pd.DataFrame, age: int,
                                           k: int = N_NEIGHBORS) -> np.ndarray:
    """Average per-nucleus kNN distances across the three two-WT-mouse references."""
    reference_rows = leave_one_wt_mouse_reference_rows(meta, age)
    distances = [mean_knn_distance(query_points, full_points[rows], k=k)
                 for _heldout, rows in reference_rows]
    return np.mean(np.stack(distances, axis=0), axis=0)


def leave_one_wt_mouse_reference_centers(full_points: np.ndarray, meta: pd.DataFrame,
                                         age: int) -> np.ndarray:
    """Make the three same-age centres, each from exactly two WT mice."""
    centers_by_mouse = leave_one_wt_mouse_reference_centers_by_mouse(full_points, meta, age)
    return np.asarray([centers_by_mouse[mouse] for mouse in sorted(centers_by_mouse)], dtype=np.float64)


def leave_one_wt_mouse_reference_centers_by_mouse(full_points: np.ndarray, meta: pd.DataFrame,
                                                  age: int) -> dict[str, np.ndarray]:
    """Return each WT mouse's age-specific centre calculated without that mouse."""
    centers = {}
    for _heldout, rows in leave_one_wt_mouse_reference_rows(meta, age):
        heldout = _heldout
        per_mouse_medians, _ = _mouse_medians(full_points[rows], meta.iloc[rows].reset_index(drop=True))
        assert len(per_mouse_medians) == 2
        centers[heldout] = per_mouse_medians.mean(axis=0)
    return centers


def shift_to_equal_loo_centers(query_points: np.ndarray, loo_centers: np.ndarray) -> float:
    query_median = np.median(np.asarray(query_points, dtype=np.float64), axis=0)
    centers = np.asarray(loo_centers, dtype=np.float64)
    assert centers.shape[0] == 3 and centers.shape[1] == query_median.shape[0]
    return float(np.mean(np.linalg.norm(centers - query_median[None, :], axis=1)))


def wt_discovery_crossfit_mask(meta: pd.DataFrame) -> np.ndarray:
    return meta.genotype.eq('WT').to_numpy() & meta.age_weeks.astype(int).isin(DISCOVERY_AGES).to_numpy()


def _choose_full_or_crossfit(full_values: np.ndarray, crossfit_values: np.ndarray,
                             use_crossfit: np.ndarray) -> np.ndarray:
    full = np.asarray(full_values, dtype=np.float64)
    crossfit = np.asarray(crossfit_values, dtype=np.float64)
    mask = np.asarray(use_crossfit, dtype=bool)
    assert full.shape == crossfit.shape and full.shape[0] == len(mask)
    selected = full.copy()
    selected[mask] = crossfit[mask]
    return selected


def assert_primary_readout_basis(mouse: pd.DataFrame, nucleus: Optional[pd.DataFrame] = None) -> None:
    """Fail if a primary WT-discovery readout was replaced by its full-fit value."""
    mouse_mask = wt_discovery_crossfit_mask(mouse)
    mouse_pairs = {
        'median_position_t': ('median_position_t_full_fit', 'median_position_t_crossfit'),
        'median_d': ('median_d_full_fit', 'median_d_crossfit'),
        'median_DPT': ('median_DPT_full_fit', 'median_DPT_crossfit'),
        'shift': ('shift_full_fit', 'shift_crossfit'),
        'spread': ('spread_full_fit', 'spread_crossfit'),
        'outside_share': ('outside_share_full_fit', 'outside_share_crossfit'),
    }
    for primary, (full_name, crossfit_name) in mouse_pairs.items():
        expected = _choose_full_or_crossfit(mouse[full_name].to_numpy(float),
                                            mouse[crossfit_name].to_numpy(float), mouse_mask)
        assert np.allclose(mouse[primary].to_numpy(float), expected, rtol=1e-12, atol=1e-12,
                           equal_nan=True), f'{primary} does not use the required WT cross-fit / KI and 12 W full-fit basis'
    if nucleus is None:
        return
    nucleus_mask = wt_discovery_crossfit_mask(nucleus)
    nucleus_pairs = {
        'position_t': ('position_t_full_fit', 'position_t_crossfit'),
        'distance_d': ('distance_d_full_fit', 'distance_d_crossfit'),
        'outside_mean_knn_distance': ('outside_mean_knn_distance_full_fit', 'outside_mean_knn_distance_crossfit'),
        'DPT': ('DPT_full_fit', 'DPT_crossfit'),
    }
    for primary, (full_name, crossfit_name) in nucleus_pairs.items():
        expected = _choose_full_or_crossfit(nucleus[full_name].to_numpy(float),
                                            nucleus[crossfit_name].to_numpy(float), nucleus_mask)
        assert np.allclose(nucleus[primary].to_numpy(float), expected, rtol=1e-12, atol=1e-12,
                           equal_nan=True), f'nucleus {primary} does not use the required fit basis'


def add_wt_reference_offsets(mouse: pd.DataFrame) -> pd.DataFrame:
    result = mouse.copy()
    result['position_offset'] = np.nan
    result['DPT_offset'] = np.nan
    for age in READOUT_AGES:
        group = result[result.age_weeks == age]
        wt = group[group.genotype == 'WT']
        if len(wt):
            rows = result.age_weeks.eq(age)
            result.loc[rows, 'position_offset'] = result.loc[rows, 'median_position_t'] - float(wt.median_position_t.mean())
            result.loc[rows, 'DPT_offset'] = result.loc[rows, 'median_DPT'] - float(wt.median_DPT.mean())
    return result


def _finite_median(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    return float(np.median(finite)) if len(finite) else np.nan


def _outside_share(values: np.ndarray, threshold: float) -> float:
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    return float(np.mean(finite > threshold)) if len(finite) and np.isfinite(threshold) else np.nan


def _weighted_metric_values(crossfit_t: dict[str, np.ndarray], crossfit_d: dict[str, np.ndarray],
                            crossfit_outside: dict[str, np.ndarray], crossfit_dpt: dict[str, np.ndarray],
                            crossfit_coords: dict[str, np.ndarray], crossfit_centers: dict[str, np.ndarray],
                            full_t: np.ndarray, full_d: np.ndarray, full_coords: np.ndarray,
                            full_dpt: np.ndarray, all_meta: pd.DataFrame,
                            outside_cutoffs: dict[int, float], full_age_centers: np.ndarray,
                            age12_center: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    medians: list[dict] = []
    per_cell: list[dict] = []
    n = len(all_meta)
    mouse_ids = all_meta.mouse_id.astype(str).to_numpy()
    ages = all_meta.age_weeks.to_numpy(dtype=int)
    genotypes = all_meta.genotype.to_numpy()
    use_crossfit = wt_discovery_crossfit_mask(all_meta)
    full_t = np.asarray(full_t, dtype=np.float64)
    full_d = np.asarray(full_d, dtype=np.float64)
    full_coords = np.asarray(full_coords, dtype=np.float64)
    full_dpt = np.asarray(full_dpt, dtype=np.float64)
    t_cf = np.full(n, np.nan); d_cf = np.full(n, np.nan); dpt_cf = np.full(n, np.nan)
    coords_cf = np.full_like(full_coords, np.nan)
    outside_cf = np.full(n, np.nan)
    for mouse, values in crossfit_t.items():
        rows = np.flatnonzero(mouse_ids == mouse)
        t_cf[rows] = values
        d_cf[rows] = crossfit_d[mouse]
        dpt_cf[rows] = crossfit_dpt[mouse]
        coords_cf[rows] = crossfit_coords[mouse]
        outside_cf[rows] = crossfit_outside[mouse]
    assert np.array_equal(np.flatnonzero(use_crossfit), np.flatnonzero(~np.isnan(t_cf)))
    assert not np.isnan(d_cf[use_crossfit]).any() and not np.isnan(dpt_cf[use_crossfit]).any()
    assert not np.isnan(outside_cf[use_crossfit]).any() and np.isfinite(coords_cf[use_crossfit]).all()
    t_values = _choose_full_or_crossfit(full_t, t_cf, use_crossfit)
    d_values = _choose_full_or_crossfit(full_d, d_cf, use_crossfit)
    dpt_values = _choose_full_or_crossfit(full_dpt, dpt_cf, use_crossfit)
    coords = _choose_full_or_crossfit(full_coords, coords_cf, use_crossfit)

    # Primary WT discovery distances use the held-mouse fold. KI distances are
    # averaged over the three same-age references, each containing two WT mice.
    outside_full = np.full(n, np.nan)
    outside_three_wt = np.full(n, np.nan)
    for age in DISCOVERY_AGES:
        ref_all = np.flatnonzero((ages == age) & (genotypes == 'WT'))
        for mouse in sorted(set(mouse_ids[ref_all])):
            held = np.flatnonzero((ages == age) & (genotypes == 'WT') & (mouse_ids == mouse))
            ref = np.flatnonzero((ages == age) & (genotypes == 'WT') & (mouse_ids != mouse))
            outside_full[held] = mean_knn_distance(full_coords[held], full_coords[ref])
        targets = np.flatnonzero((ages == age) & (genotypes == 'KI'))
        if len(targets):
            outside_full[targets] = mean_knn_distance_equal_loo_references(
                full_coords[targets], full_coords, all_meta, age)
        all_age_rows = np.flatnonzero(ages == age)
        outside_three_wt[all_age_rows] = mean_knn_distance(full_coords[all_age_rows], full_coords[ref_all])
    age12_wt = np.flatnonzero((ages == 12) & (genotypes == 'WT'))
    age12_ki = np.flatnonzero((ages == 12) & (genotypes == 'KI'))
    for mouse in sorted(set(mouse_ids[age12_wt])):
        held = np.flatnonzero((ages == 12) & (genotypes == 'WT') & (mouse_ids == mouse))
        ref = np.flatnonzero((ages == 12) & (genotypes == 'WT') & (mouse_ids != mouse))
        outside_full[held] = mean_knn_distance(full_coords[held], full_coords[ref])
    if len(age12_ki):
        outside_full[age12_ki] = mean_knn_distance_equal_loo_references(
            full_coords[age12_ki], full_coords, all_meta, 12)
    age12_rows = np.flatnonzero(ages == 12)
    outside_three_wt[age12_rows] = mean_knn_distance(full_coords[age12_rows], full_coords[age12_wt])
    outside_values = _choose_full_or_crossfit(outside_full, outside_cf, use_crossfit)
    loo_wt_centers_by_mouse = {
        age: leave_one_wt_mouse_reference_centers_by_mouse(full_coords, all_meta, age)
        for age in READOUT_AGES
    }

    for mouse in sorted(set(mouse_ids)):
        rows = np.flatnonzero(mouse_ids == mouse)
        first = all_meta.iloc[rows[0]]
        age = int(first.age_weeks)
        is_crossfit = bool(use_crossfit[rows].all())
        assert bool(use_crossfit[rows].any()) == is_crossfit, (mouse, 'mixed per-nucleus fit basis within a mouse')
        if age == 12:
            full_center = age12_center
        else:
            full_center = full_age_centers[DISCOVERY_AGES.index(age)]
        cf_center = crossfit_centers.get(mouse)
        full_median_coord = np.median(full_coords[rows], axis=0)
        full_three_wt_center_shift = float(np.linalg.norm(full_median_coord - full_center))
        if first.genotype == 'KI':
            full_shift_equal_reference = shift_to_equal_loo_centers(
                full_coords[rows], np.asarray(list(loo_wt_centers_by_mouse[age].values())))
        elif age == 12:
            # Match each 12 W WT mouse to the same two-mouse leave-one-out
            # reference used for its outside distance; never include itself.
            full_shift_equal_reference = float(np.linalg.norm(
                full_median_coord - loo_wt_centers_by_mouse[age][mouse]))
        else:
            full_shift_equal_reference = full_three_wt_center_shift
        full_shift = full_shift_equal_reference
        full_spread = float(np.median(np.linalg.norm(full_coords[rows] - full_median_coord, axis=1)))
        if is_crossfit:
            cf_median_coord = np.median(coords_cf[rows], axis=0)
            cf_shift = float(np.linalg.norm(cf_median_coord - cf_center))
            cf_spread = float(np.median(np.linalg.norm(coords_cf[rows] - cf_median_coord, axis=1)))
        else:
            cf_shift = np.nan
            cf_spread = np.nan
        threshold = outside_cutoffs.get(age, np.nan)
        full_outside_share = _outside_share(outside_full[rows], threshold)
        cf_outside_share = _outside_share(outside_cf[rows], threshold) if is_crossfit else np.nan
        full_dpt_median = _finite_median(full_dpt[rows])
        cf_dpt_median = _finite_median(dpt_cf[rows]) if is_crossfit else np.nan
        readout_basis = ('WT_discovery_leave_one_mouse_out' if is_crossfit else
                         ('full_WT_discovery_fit_equal_two_WT_references' if first.genotype == 'KI'
                          else 'full_WT_discovery_fit'))
        record = dict(
            mouse_id=mouse, age_weeks=age, genotype=str(first.genotype), sex=str(first.sex), n_nuclei=len(rows),
            readout_basis=readout_basis,
            median_position_t=_finite_median(t_values[rows]),
            median_position_t_full_fit=_finite_median(full_t[rows]),
            median_position_t_crossfit=_finite_median(t_cf[rows]) if is_crossfit else np.nan,
            median_d=_finite_median(d_values[rows]), median_d_full_fit=_finite_median(full_d[rows]),
            median_d_crossfit=_finite_median(d_cf[rows]) if is_crossfit else np.nan,
            median_DPT=_finite_median(dpt_values[rows]), median_DPT_full_fit=full_dpt_median,
            median_DPT_crossfit=cf_dpt_median,
            DPT_infinite_n=int(np.isinf(dpt_values[rows]).sum()),
            DPT_full_fit_infinite_n=int(np.isinf(full_dpt[rows]).sum()),
            DPT_crossfit_infinite_n=int(np.isinf(dpt_cf[rows]).sum()) if is_crossfit else 0,
            shift=cf_shift if is_crossfit else full_shift,
            shift_full_fit=full_shift, shift_crossfit=cf_shift,
            shift_full_three_WT_center=full_three_wt_center_shift,
            spread=float(np.median(np.linalg.norm(coords[rows] - np.median(coords[rows], axis=0), axis=1))),
            spread_full_fit=full_spread, spread_crossfit=cf_spread,
            outside_share=_outside_share(outside_values[rows], threshold),
            outside_share_full_fit=full_outside_share, outside_share_crossfit=cf_outside_share,
            position_12w_between_5_and_18=bool(5 <= np.median(t_values[rows]) <= 18)
            if age == 12 and first.genotype == 'WT' else np.nan)
        medians.append(record)
        for i in rows:
            per_cell.append(dict(
                nucleus_id=str(all_meta.iloc[i].nucleus_id), mouse_id=mouse, age_weeks=age,
                genotype=str(first.genotype), sex=str(first.sex),
                readout_basis=readout_basis,
                position_t=float(t_values[i]), position_t_full_fit=float(full_t[i]),
                position_t_crossfit=float(t_cf[i]), distance_d=float(d_values[i]),
                distance_d_full_fit=float(full_d[i]), distance_d_crossfit=float(d_cf[i]),
                outside_mean_knn_distance=float(outside_values[i]),
                outside_mean_knn_distance_full_fit=float(outside_full[i]),
                outside_mean_knn_distance_crossfit=float(outside_cf[i]),
                outside_mean_knn_distance_three_WT_reference=float(outside_three_wt[i]),
                outside_threshold_p95=float(threshold),
                outside_flag=bool(outside_values[i] > threshold)
                if np.isfinite(outside_values[i]) and np.isfinite(threshold) else np.nan,
                DPT=float(dpt_values[i]), DPT_full_fit=float(full_dpt[i]),
                DPT_crossfit=float(dpt_cf[i]), DPT_infinite=bool(np.isinf(dpt_values[i])),
                DPT_full_fit_infinite=bool(np.isinf(full_dpt[i])),
                DPT_crossfit_infinite=bool(np.isinf(dpt_cf[i])) if is_crossfit else np.nan))
    mouse = pd.DataFrame(medians)
    nucleus = pd.DataFrame(per_cell)
    assert_primary_readout_basis(mouse, nucleus)
    return mouse, nucleus


def exact_two_sided(values_ki: np.ndarray, values_wt: np.ndarray) -> tuple[float, float, str]:
    ki = np.asarray(values_ki, dtype=np.float64); wt = np.asarray(values_wt, dtype=np.float64)
    if len(ki) != 3 or len(wt) != 3 or not np.isfinite(ki).all() or not np.isfinite(wt).all():
        return np.nan, np.nan, 'NOT_EVALUABLE'
    observed = float(ki.mean() - wt.mean())
    all_values = np.concatenate((ki, wt))
    diffs = []
    for pick in itertools.combinations(range(6), 3):
        is_ki = np.zeros(6, dtype=bool); is_ki[list(pick)] = True
        diffs.append(float(all_values[is_ki].mean() - all_values[~is_ki].mean()))
    p = sum(abs(v) >= abs(observed) - 1e-12 for v in diffs) / 20.0
    if np.min(ki) > np.max(wt): sep = 'KI_GT_WT'
    elif np.max(ki) < np.min(wt): sep = 'KI_LT_WT'
    else: sep = 'none'
    return observed, float(p), sep


def make_age_tests(mouse: pd.DataFrame, metrics: tuple[str, ...]) -> tuple[pd.DataFrame, pd.DataFrame]:
    age_rows = []
    flags = []
    for metric in metrics:
        for age in READOUT_AGES:
            group = mouse[mouse.age_weeks == age]
            ki = group[group.genotype == 'KI'][metric].to_numpy(dtype=np.float64)
            wt = group[group.genotype == 'WT'][metric].to_numpy(dtype=np.float64)
            difference, p, separation = exact_two_sided(ki, wt)
            age_rows.append(dict(age_weeks=age, measure=metric, n_KI=len(ki), n_WT=len(wt),
                                 KI_minus_WT_mean=difference, complete_separation=separation,
                                 exact_two_sided_permutation_p=p, permutation_labellings=20))
        discovery = [r for r in age_rows if r['measure'] == metric and r['age_weeks'] in DISCOVERY_AGES]
        signs = [int(np.sign(r['KI_minus_WT_mean'])) if np.isfinite(r['KI_minus_WT_mean']) else 0 for r in discovery]
        pos = sum(s > 0 for s in signs); neg = sum(s < 0 for s in signs)
        direction = 1 if pos > neg else (-1 if neg > pos else 0)
        sep_label = 'KI_GT_WT' if direction > 0 else ('KI_LT_WT' if direction < 0 else '')
        same_sign_n = max(pos, neg)
        sep_n = sum(r['complete_separation'] == sep_label for r in discovery) if sep_label else 0
        replicated_sign = 0
        age12 = next(r for r in age_rows if r['measure'] == metric and r['age_weeks'] == 12)
        if np.isfinite(age12['KI_minus_WT_mean']):
            replicated_sign = int(np.sign(age12['KI_minus_WT_mean']))
        flags.append(dict(measure=metric, discovery_direction=('KI_GT_WT' if direction > 0 else ('KI_LT_WT' if direction < 0 else 'none')),
                          discovery_same_sign_ages=same_sign_n, discovery_complete_separation_same_direction_ages=sep_n,
                          consistent=bool(direction != 0 and same_sign_n >= 3 and sep_n >= 2),
                          replicated_12w_same_sign=bool(direction != 0 and replicated_sign == direction),
                          age12_direction=('KI_GT_WT' if replicated_sign > 0 else ('KI_LT_WT' if replicated_sign < 0 else 'tie_or_not_evaluable'))))
    return pd.DataFrame(age_rows), pd.DataFrame(flags)


def _fit_one_scope(X: np.ndarray, meta: pd.DataFrame, train_mask: np.ndarray, query_rows: np.ndarray,
                   seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    model = fit_wt_model(X, meta, train_mask, n_pcs=N_PCS, seed=seed)
    train_rows = np.flatnonzero(train_mask)
    train_pc = transform_pcs(X, train_rows, model)
    train_meta = meta.iloc[train_rows].reset_index(drop=True)
    # Path and DPT operate in the fitted PCA score coordinates; no whitening or
    # post-PCA rescaling is part of the frozen representation contract.
    train_path_points = train_pc
    centers = age_centers(train_path_points, train_meta)
    query_pc = transform_pcs(X, query_rows, model)
    query_path_points = query_pc
    t, d, foot = project_to_path(query_path_points, centers)
    return t, d, foot, query_path_points, dict(model=model, centers=centers,
                                               train_pc=train_pc, train_meta=train_meta,
                                               train_path_points=train_path_points,
                                               query_pc=query_pc, query_path_points=query_path_points)


def analyze_space(X: np.ndarray, meta: pd.DataFrame, out_dir: Path, space: str, population: str,
                  pca_diagnostic_n: int = PCA_DIAGNOSTIC_NUCLEI) -> dict:
    """Fit WT discovery representation, cross-fit WT discovery mice, and compute readouts."""
    meta = meta.reset_index(drop=True).copy()
    meta['mouse_id'] = meta.mouse_id.astype(str)
    assert X.dtype == np.float64 and X.ndim == 2 and X.shape[0] == len(meta)
    assert np.isfinite(meta.log10_total_umi.to_numpy(dtype=float)).all()
    assert (meta.n_detected.to_numpy(dtype=int) >= 0).all()
    assert meta.genotype.isin(['WT', 'KI']).all()
    assert set(meta.age_weeks.astype(int)) == set(READOUT_AGES)
    mouse_group_sizes = meta.groupby('mouse_id')[['age_weeks', 'genotype', 'sex']].nunique()
    assert mouse_group_sizes.le(1).all().all(), 'a mouse ID maps to multiple ages, genotypes, or sexes'
    for age in READOUT_AGES:
        for genotype in ('WT', 'KI'):
            n_mice = meta.loc[(meta.age_weeks == age) & (meta.genotype == genotype), 'mouse_id'].nunique()
            assert n_mice == 3, (age, genotype, n_mice, 'expected exactly three mice')
    fit_mask = meta.genotype.eq('WT').to_numpy() & meta.age_weeks.isin(DISCOVERY_AGES).to_numpy()
    fit_meta = meta.loc[fit_mask].reset_index(drop=True)
    counts = fit_meta.groupby(['age_weeks', 'mouse_id']).size().groupby(level=0).size()
    assert counts.to_dict() == {age: 3 for age in DISCOVERY_AGES}, ('expected three WT discovery mice at each age', counts.to_dict())
    assert fit_meta.mouse_id.nunique() == 12
    all_rows = np.arange(len(meta), dtype=int)
    full_t, full_d, _full_foot, full_coords, fitted = _fit_one_scope(X, meta, fit_mask, all_rows, DPT_SEED)
    full_centers = fitted['centers']
    pca_diag_table, pca_diag_summary = pca_diagnostics(
        X, meta, fit_mask, fitted['model'], diagnostic_n=pca_diagnostic_n)
    full_dpt, full_dpt_info = dpt_from_wt_graph(fitted['train_path_points'], fitted['train_meta'],
                                                full_coords, meta, seed=DPT_SEED)
    graph_node_rows = []
    for node in full_dpt_info.pop('graph_nodes'):
        node.update(space=space, population=population, scope='full_WT_discovery_graph',
                    heldout_mouse='NONE', root_nucleus_id=full_dpt_info['root_nucleus_id'])
        graph_node_rows.append(node)

    cross_t: dict[str, np.ndarray] = {}
    cross_d: dict[str, np.ndarray] = {}
    cross_out: dict[str, np.ndarray] = {}
    cross_dpt: dict[str, np.ndarray] = {}
    cross_coords: dict[str, np.ndarray] = {}
    cross_centers: dict[str, np.ndarray] = {}
    inf_rows = []
    for mouse_id in sorted(fit_meta.mouse_id.unique()):
        held_rows = np.flatnonzero(fit_mask & meta.mouse_id.eq(mouse_id).to_numpy())
        train_mask = fit_mask & ~meta.mouse_id.eq(mouse_id).to_numpy()
        t, d, _foot, held_coords, fold = _fit_one_scope(X, meta, train_mask, held_rows, DPT_SEED)
        held_meta = meta.iloc[held_rows].reset_index(drop=True)
        age = int(held_meta.age_weeks.iloc[0])
        cross_t[mouse_id] = t
        cross_d[mouse_id] = d
        cross_coords[mouse_id] = held_coords
        cross_centers[mouse_id] = fold['centers'][DISCOVERY_AGES.index(age)]
        train_same_age = fold['train_meta'].age_weeks.to_numpy(dtype=int) == age
        cross_out[mouse_id] = mean_knn_distance(held_coords, fold['train_path_points'][train_same_age])
        loo_dpt, loo_info = dpt_from_wt_graph(fold['train_path_points'], fold['train_meta'],
                                               fold['query_path_points'], held_meta, seed=DPT_SEED)
        cross_dpt[mouse_id] = loo_dpt
        for node in loo_info.pop('graph_nodes'):
            node.update(space=space, population=population, scope='WT_discovery_leave_one_mouse_out_graph',
                        heldout_mouse=mouse_id, root_nucleus_id=loo_info['root_nucleus_id'])
            graph_node_rows.append(node)
        inf_rows.append(dict(scope='WT_discovery_leave_one_mouse_out', mouse_id=mouse_id,
                             n_infinite=int(np.isinf(loo_dpt).sum()), n_nuclei=len(loo_dpt),
                             graph_infinite=int(loo_info['n_graph_infinite']), n_graph_nodes=loo_info['n_graph_nodes'],
                             root_nucleus_id=loo_info['root_nucleus_id']))

    # Age 12 is held out from both the path and the WT discovery graph. Its WT outside-distance
    # calibration is leave-one-mouse-out among the three WT 12 W mice.
    ages = meta.age_weeks.to_numpy(dtype=int)
    mice = meta.mouse_id.astype(str).to_numpy()
    genotype = meta.genotype.to_numpy()
    age12_wt_rows = np.flatnonzero((ages == 12) & (genotype == 'WT'))
    age12_cross_values = []
    for mouse_id in sorted(set(mice[age12_wt_rows])):
        held = np.flatnonzero((ages == 12) & (genotype == 'WT') & (mice == mouse_id))
        reference = np.flatnonzero((ages == 12) & (genotype == 'WT') & (mice != mouse_id))
        values = mean_knn_distance(full_coords[held], full_coords[reference])
        age12_cross_values.extend(values.tolist())

    outside_cutoffs = {}
    for age in DISCOVERY_AGES:
        age_mice = fit_meta.loc[fit_meta.age_weeks.eq(age), 'mouse_id'].astype(str).unique()
        values = np.concatenate([cross_out[mouse_id] for mouse_id in sorted(age_mice)])
        outside_cutoffs[age] = float(np.quantile(values, .95, method='linear'))
    outside_cutoffs[12] = float(np.quantile(np.asarray(age12_cross_values), .95, method='linear'))
    age12_medians, _ = _mouse_medians(full_coords[age12_wt_rows], meta.iloc[age12_wt_rows].reset_index(drop=True))
    age12_center = age12_medians.mean(axis=0)
    mouse, nucleus = _weighted_metric_values(
        cross_t, cross_d, cross_out, cross_dpt, cross_coords, cross_centers,
        full_t, full_d, full_coords, full_dpt, meta, outside_cutoffs, full_centers, age12_center)

    discovery_mouse = mouse[(mouse.age_weeks.isin(DISCOVERY_AGES)) & (mouse.genotype == 'WT')]
    dpt_rho_result = spearmanr(discovery_mouse.age_weeks.to_numpy(dtype=float),
                               discovery_mouse.median_DPT_crossfit.to_numpy(dtype=float), nan_policy='omit')
    dpt_rho = float(dpt_rho_result.statistic) if np.isfinite(dpt_rho_result.statistic) else np.nan
    pos_rho_result = spearmanr(discovery_mouse.age_weeks.to_numpy(dtype=float),
                               discovery_mouse.median_position_t_crossfit.to_numpy(dtype=float), nan_policy='omit')
    pos_rho = float(pos_rho_result.statistic) if np.isfinite(pos_rho_result.statistic) else np.nan

    mouse = add_wt_reference_offsets(mouse)
    assert_primary_readout_basis(mouse, nucleus)
    metrics = ('position_offset', 'median_d', 'shift', 'spread', 'outside_share', 'DPT_offset')
    age_tests, reading_flags = make_age_tests(mouse, metrics)
    pos_valid = bool(np.isfinite(pos_rho) and pos_rho >= SPEARMAN_FLOOR)
    dpt_valid = bool(np.isfinite(dpt_rho) and dpt_rho >= SPEARMAN_FLOOR)
    mouse['position_informative'] = pos_valid
    mouse['DPT_informative'] = dpt_valid
    validity = pd.DataFrame([
        dict(space=space, population=population, readout='position_t',
             spearman_rho_age_vs_crossfit_WT_mouse_median=pos_rho,
             threshold=SPEARMAN_FLOOR, informative=pos_valid,
             WT_12w_mice_between_5_and_18=int(mouse.loc[(mouse.age_weeks == 12) & (mouse.genotype == 'WT'), 'position_12w_between_5_and_18'].sum())),
        dict(space=space, population=population, readout='DPT',
             spearman_rho_age_vs_crossfit_WT_mouse_median=dpt_rho,
             threshold=SPEARMAN_FLOOR, informative=dpt_valid, WT_12w_mice_between_5_and_18=np.nan)])
    for row in inf_rows:
        row.update(space=space, population=population, scanpy_version=full_dpt_info['scanpy_version'])
    for mouse_id in sorted(meta.mouse_id.astype(str).unique()):
        rows = np.flatnonzero(meta.mouse_id.astype(str).to_numpy() == mouse_id)
        inf_rows.append(dict(space=space, population=population,
                             scope='full_WT_discovery_graph_query_by_mouse', mouse_id=mouse_id,
                             n_infinite=int(np.isinf(full_dpt[rows]).sum()), n_nuclei=len(rows),
                             graph_infinite=np.nan, n_graph_nodes=np.nan,
                             root_nucleus_id=full_dpt_info['root_nucleus_id'],
                             scanpy_version=full_dpt_info['scanpy_version']))
    inf_rows.append(dict(space=space, population=population,
                         scope='full_WT_discovery_graph_and_all_transfers', mouse_id='ALL',
                         n_infinite=int(np.isinf(full_dpt).sum()), n_nuclei=len(full_dpt),
                         graph_infinite=int(full_dpt_info['n_graph_infinite']),
                         n_graph_nodes=full_dpt_info['n_graph_nodes'], root_nucleus_id=full_dpt_info['root_nucleus_id'],
                         scanpy_version=full_dpt_info['scanpy_version']))

    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f'{space}_{POP_SLUG[population]}'
    mouse.insert(0, 'population', population); mouse.insert(0, 'space', space)
    nucleus.insert(0, 'population', population); nucleus.insert(0, 'space', space)
    mouse.to_csv(out_dir / f'{tag}_MOUSE_READOUT.tsv', sep='\t', index=False)
    nucleus.to_csv(out_dir / f'{tag}_NUCLEUS_READOUT.tsv.gz', sep='\t', index=False, compression='gzip')
    age_tests.insert(0, 'population', population); age_tests.insert(0, 'space', space)
    reading_flags.insert(0, 'population', population); reading_flags.insert(0, 'space', space)
    age_tests.to_csv(out_dir / f'{tag}_AGE_TESTS.tsv', sep='\t', index=False)
    reading_flags.to_csv(out_dir / f'{tag}_READING_FLAGS.tsv', sep='\t', index=False)
    validity.to_csv(out_dir / f'{tag}_VALIDITY.tsv', sep='\t', index=False)
    pd.DataFrame(inf_rows).to_csv(out_dir / f'{tag}_DPT_INF_COUNTS.tsv', sep='\t', index=False)
    pd.DataFrame(graph_node_rows).to_csv(out_dir / f'{tag}_DPT_GRAPH_NODES.tsv.gz', sep='\t', index=False,
                                         compression='gzip')
    pca_diag_table.to_csv(out_dir / f'{tag}_PCA_DIAGNOSTICS.tsv', sep='\t', index=False)
    np.savez_compressed(out_dir / f'{tag}_FULL_MODEL.npz', slope=fitted['model']['slope'],
                        intercept=fitted['model']['intercept'], residual_mean=fitted['model']['residual_mean'],
                        residual_sd=fitted['model']['residual_sd'], z_mean=fitted['model']['z_mean'],
                        rotation=fitted['model']['rotation'], eigenvalues=fitted['model']['eigenvalues'],
                        variance_fractions=fitted['model']['variance_fractions'],
                        cumulative_variance_fractions=fitted['model']['cumulative_variance_fractions'],
                        total_weighted_variance=fitted['model']['total_weighted_variance'],
                        pca_diagnostic_subsample_seed=PCA_DIAGNOSTIC_SEED,
                        age_centers=full_centers)
    return dict(space=space, population=population, n_nuclei=len(meta), n_features=int(X.shape[1]),
                n_pcs=int(fitted['model']['n_pcs']), n_WT_discovery_mice=12,
                position_spearman_rho=pos_rho, position_informative=pos_valid,
                DPT_spearman_rho=dpt_rho, DPT_informative=dpt_valid,
                pca_diagnostics=pca_diag_summary,
                scanpy_version=full_dpt_info['scanpy_version'], DPT_graph=full_dpt_info)


def _interface_sum_map(path: Path) -> dict[str, str]:
    entries = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        assert len(parts) == 2 and len(parts[0]) == 64, ('malformed SHA256SUMS entry', line)
        name = parts[1].strip().lstrip('*')
        entries[Path(name).name] = parts[0]
    return entries


def _read_interface(input_dir: Path, space: str, population: str) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame, dict]:
    tag = f'{space}_{POP_SLUG[population]}'
    matrix_path = input_dir / f'state_{tag}.npy'
    nuclei_path = input_dir / f'state_{tag}_nuclei.tsv'
    features_path = input_dir / f'state_{tag}_features.tsv'
    assert matrix_path.is_file() and nuclei_path.is_file() and features_path.is_file(), f'missing interface triplet for {tag}'
    sums_path = input_dir / 'SHA256SUMS.txt'
    assert sums_path.is_file(), f'missing checksum file for interface {input_dir}'
    listed = _interface_sum_map(sums_path)
    for source in (matrix_path, nuclei_path, features_path):
        expected = listed.get(source.name)
        assert expected is not None, f'{source.name} is absent from the interface SHA256SUMS.txt'
        actual = sha(source)
        assert actual == expected, (f'{source.name} differs from interface SHA256SUMS.txt', actual, expected)
    matrix = np.load(matrix_path, mmap_mode='r', allow_pickle=False)
    nuclei = pd.read_csv(nuclei_path, sep='\t', dtype={'nucleus_id': str, 'mouse_id': str, 'genotype': str, 'sex': str})
    features = pd.read_csv(features_path, sep='\t', dtype={'stable_id': str})
    required = {'nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex', 'log10_total_umi', 'n_detected'}
    assert required <= set(nuclei.columns), (tag, required - set(nuclei.columns))
    assert 'stable_id' in features.columns
    if 'feature_index' not in features.columns:
        features.insert(0, 'feature_index', np.arange(len(features), dtype=int))
    assert matrix.dtype == np.float64 and matrix.ndim == 2 and matrix.shape == (len(nuclei), len(features))
    assert not features.feature_index.duplicated().any()
    assert not nuclei.nucleus_id.duplicated().any() and not features.stable_id.duplicated().any()
    assert nuclei.age_weeks.isin(READOUT_AGES).all()
    assert np.isfinite(nuclei[['log10_total_umi', 'n_detected']].to_numpy(dtype=float)).all()
    for first in range(0, len(matrix), CHUNK):
        assert np.isfinite(matrix[first:first + CHUNK]).all(), f'{tag}: matrix contains nonfinite values'
    return matrix, nuclei, features, {str(p): sha(p) for p in (matrix_path, nuclei_path, features_path)}


def run_real(p_dir: Path, g_dir: Path, out_dir: Path, gate_path: Path) -> dict:
    # Gate/sums hashes are checked before any state matrix is opened.
    gate = check_go(gate_path, expected_scope='shared_readout',
                    expected_scripts={Path(__file__).name: sha(Path(__file__))},
                    p_dir=p_dir, g_dir=g_dir)
    assert not out_dir.exists(), f'output exists; never overwrite: {out_dir}'
    env_record, explicit = validate_env_provenance(capture=True)
    out_dir.mkdir(parents=True)
    (out_dir / 'CONDA_EXPLICIT.txt').write_bytes(explicit)
    started = time.time()
    records = []; input_hashes = {}; metadata_by_pop = {}; metadata_checks = []
    for space in ('P', 'G_a', 'G_b'):
        for pop in POPS:
            input_dir = p_dir if space == 'P' else g_dir
            matrix, nuclei, features, hashes = _read_interface(input_dir, space, pop)
            input_hashes.update(hashes)
            signature_columns = ['nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex',
                                 'log10_total_umi', 'n_detected']
            signature = nuclei[signature_columns].sort_values('nucleus_id').reset_index(drop=True)
            if pop not in metadata_by_pop:
                metadata_by_pop[pop] = signature
            reference = metadata_by_pop[pop]
            exact_columns = ['nucleus_id', 'mouse_id', 'age_weeks', 'genotype', 'sex', 'n_detected']
            pd.testing.assert_frame_equal(reference[exact_columns], signature[exact_columns],
                                          check_dtype=False, check_exact=True, check_like=False)
            a = reference.log10_total_umi.to_numpy(dtype=np.float64)
            b = signature.log10_total_umi.to_numpy(dtype=np.float64)
            total_a, total_b = 10 ** a, 10 ** b
            rounded_a, rounded_b = np.rint(total_a), np.rint(total_b)
            for implied, rounded in ((total_a, rounded_a), (total_b, rounded_b)):
                assert np.isfinite(implied).all() and (rounded > 0).all(), 'invalid implied UMI total'
                assert (np.abs(implied - rounded) <= 1e-6 * rounded).all(), 'implied UMI total is not an integer'
            assert np.array_equal(rounded_a, rounded_b), 'cross-space integer UMI totals differ'
            delta = np.abs(a - b)
            max_delta = float(np.max(delta))
            assert max_delta <= 1e-12, 'cross-space log10 UMI delta exceeds rounding allowance'
            metadata_checks.append(dict(population=pop, space=space, reference_space='P',
                                        n_nuclei=len(signature), n_log10_rows_differing=int(np.count_nonzero(a != b)),
                                        max_abs_log10_delta=max_delta,
                                        integer_total_relative_tolerance=1e-6, log10_absolute_tolerance=1e-12))
            records.append(analyze_space(matrix, nuclei, out_dir, space, pop))
    summary = dict(mode='real', spaces=['P', 'G_a', 'G_b'], populations=list(POPS),
                   readouts=['position_t', 'DPT'], n_scopes=len(records), scopes=records,
                    input_interface_sha256=input_hashes, environment=env_record,
                    cross_space_metadata_check=metadata_checks,
                    script_sha256={Path(__file__).name: sha(Path(__file__))},
                    go_file_sha256=sha(gate_path),
                   dpt_contract=dict(k=N_NEIGHBORS, n_dcs=N_DCS, seed=DPT_SEED, subsample='equal nuclei per WT discovery mouse, count=min mouse count',
                                     disconnected='count + exclude inf from medians', scanpy_version=ENV_VERSIONS['scanpy']),
                    elapsed_seconds=time.time() - started, source_H5AD_not_read_by_shared_module=True)
    summary = json_safe(summary)
    record = json_safe(dict(mode='real', approval_gate=gate, summary=summary,
        conda_explicit_command=env_record['export_command'], conda_explicit_sha256=env_record['explicit_sha256'],
        environment_versions=env_record['versions']))
    (out_dir / 'SUMMARY.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out_dir / 'RECORD.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    sums = []
    for p in sorted(out_dir.iterdir()):
        if p.is_file() and p.name != 'SHA256SUMS.txt': sums.append(f'{sha(p)}  {p.name}\n')
    (out_dir / 'SHA256SUMS.txt').write_text(''.join(sums), encoding='utf-8')
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--real', action='store_true', required=True)
    ap.add_argument('--p-dir', type=Path, required=True, help=' P-matrix interface directory')
    ap.add_argument('--g-dir', type=Path, required=True, help='all-gene G(a)/G(b) interface directory')
    ap.add_argument('--go-file', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=REAL_OUT)
    args = ap.parse_args()
    if not sys.flags.isolated or not sys.dont_write_bytecode:
        raise SystemExit('refusing: launch this entry point with python -I -B')
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('refusing: checks require Python assertions')
    # Scanpy/numba need the single_cell environment's own libstdc++ at process start.
    envlib = str(ENV / 'lib')
    if envlib not in os.environ.get('LD_LIBRARY_PATH', '').split(':'):
        raise SystemExit('refusing: launch with LD_LIBRARY_PATH=' + envlib + ':$LD_LIBRARY_PATH')
    try:
        result = run_real(args.p_dir, args.g_dir, args.out, args.go_file)
        print(json.dumps(dict(completed=True, scopes=len(result['scopes']), output=str(args.out))))
    except BaseException:
        if args.out.exists():
            (args.out / 'FAILED.json').write_text(json.dumps(dict(traceback=traceback.format_exc()), indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()

#%% debug
