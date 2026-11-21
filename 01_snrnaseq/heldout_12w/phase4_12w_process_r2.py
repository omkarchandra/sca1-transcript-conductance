
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# r2 (pre-specified): parity on off_manifold and state_id against the saved S07 per-nucleus table, row mapping through nucleus_metadata.tsv.gz, replay also compares the projection, FAILED.json on failure.
# S07 definition of off_manifold, STRICT '>': scripts/s07_projection.py line 40 'off_manifold=minimum>off_manifold_threshold' and scripts/s07_assemble_states.py line 122 "p['off_manifold']=p['distance']>threshold".
# Phase 4: process the six excluded 12-week mice (Purkinje cell and Bergmann glia only) EXACTLY like the 24 mice were processed in S06, then PLACE them on the FROZEN S07 discovery models. No refit of anything.
# the smoke (one mouse) reads real expression, so it is part of the run. Run with: python -I -B 01_snrnaseq/heldout_12w/phase4_12w_process_r2.py --smoke | --full | --check | --replay DIR
# Reuses the S06/S07 code instead of copying it: s06_features.normalize, p00.counts_scan.plan_chunks, s07_assemble_states.projection / assign_states. Label-blind until the very last step: genotype and sex are joined by animal_id AFTER placement.
# review's rule (): the 24 mice are in-sample for the discovery model, the 12 W mice are not. PC-score means per mouse are the PRIMARY output; distance / off_manifold / state are written but labelled out-of-sample (locked threshold, not comparable with 24-mouse values).
import argparse, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, h5py
from anndata.io import read_elem
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p00.counts_scan import plan_chunks
from s06_features import normalize
from s07_assemble_states import projection, assign_states

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
POPS = ('Purkinje cell', 'Bergmann glia')
SIX = ('2061', '2062', '22018', '2063', '2065', '22019')            # WT 2061, 2062, 22018; KI 2063, 2065, 22019
CFG = dict(h5ad=ROOT / _release_paths.expand('${DATA_ROOT}/data/adata_leiden_0.175_annotated_cell_types_v3.h5ad'), manifest=ROOT / _release_paths.expand('${DATA_ROOT}/results/input_sca1_audit_run01/nucleus_manifest.tsv.gz'),
           animal_manifest=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'), A=ROOT / _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01'), CORE=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01'),
           TREES=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_fullprecision_core_trees_v01'), S07PROJ=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02'), REGISTRY=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_state_registry_v01/STATE_REGISTRY.json'), DATE='20261007')
PINS = {  # relative path -> full sha256 (from ${DATA_ROOT}/phase4_12w_scoping/SCOPING.json, accepted by review at )
    _release_paths.expand('${DATA_ROOT}/data/adata_leiden_0.175_annotated_cell_types_v3.h5ad'): '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86',
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07',
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/feature_registry.tsv'): 'd49b138d3a806db64fb37b6994a30034569cd530bfe3b6a82becf10f367f8c24',
    _release_paths.expand('${DATA_ROOT}/s12_12w_sex_lookup/SEX_12W.json'): '0e2715e8b7ef3159c79525267734e628738d93b74e04fbbf9ef35d98839fdeb1',
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Purkinje_cell_projection.tsv.gz'): '8605f82700db2e783032ce66dd744be479af987186595a0d6d09d4a860c6b458',
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Bergmann_glia_projection.tsv.gz'): '29b779ee72e327e066e5e9efa8ece185560498b1dbca650cd71aa9019342ae10',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): '23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540',
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/cell_metadata.tsv.gz'): 'ee6c4bca63f237825d7d80389e759263dbb45477e8535abb015b350f0595dbe6'}
SCOPING_SHA256 = 'a7131e744951ec404f15ba003542a80ad793df814ea17c58f6290a05b8ec7447'   # the scoping record carries the hashes of the nucleus manifest and of every frozen S06/S07 file used below


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def verify_pins():
    """Precondition: every input and every frozen model file equals its pinned hash. Returns {path: sha256}."""
    assert sha(ROOT / _release_paths.expand('${DATA_ROOT}/phase4_12w_scoping/SCOPING.json')) == SCOPING_SHA256
    sc = json.loads((ROOT / _release_paths.expand('${DATA_ROOT}/phase4_12w_scoping/SCOPING.json')).read_text()); pins = dict(PINS)
    pins[_release_paths.expand('${DATA_ROOT}/results/input_sca1_audit_run01/nucleus_manifest.tsv.gz')] = sc['nucleus_manifest']['sha256']
    for pop in ('Purkinje_cell', 'Bergmann_glia'):
        for part in sc['frozen_S07_files'][pop].values(): pins[part['path']] = part['sha256']
    for k in ('state_registry', 'stageA_transform_keys', 'stageA_transform_parameters'): pins[sc['frozen_S07_files'][k]['path']] = sc['frozen_S07_files'][k]['sha256']
    return check_pins(pins, ROOT)


def check_pins(pins, root):
    got = {}
    for rel, h in pins.items():
        got[rel] = sha(Path(root) / rel); assert got[rel] == h, ('pinned file changed', rel)
    return got


def load_model(pop, cfg=CFG):
    name = 'discovery__' + pop.replace(' ', '_'); A = cfg['A']
    keys = json.loads((A / 'transform_keys.json').read_text()); k = next(i for i, key in enumerate(keys) if key['mode'] == 'discovery' and key['population'] == pop)
    params = np.load(A / 'transform_parameters.npz'); eligible = params['eligible'][k]
    mean, scale = params['mean'][k, eligible], params['scale'][k, eligible]
    pca = np.load(cfg['CORE'] / name / 'pca_model.npz')
    assert np.array_equal(pca['feature_indices'], np.flatnonzero(eligible)) and np.allclose(pca['mean'], mean) and np.allclose(pca['scale'], scale)
    reg = json.loads(cfg['REGISTRY'].read_text()); entry = next(s for s in reg['populations'] if s['population'] == pop)
    return dict(name=name, eligible=eligible, mean=mean, scale=scale, rotation=pca['rotation'], tree=json.loads((cfg['TREES'] / name / 'primary_tree.json').read_text()),
                entries=[s for s in reg['states'] if s['population'] == pop], bins=entry['bins'], threshold=float(entry['training_distance_threshold']))


def apply_frozen_transform(normalized_rows, model):
    _release_paths.expand('${DATA_ROOT}/helpers/s07_prepare_core_pca.py line 55, verbatim formula, with the frozen eligible/mean/scale/rotation.')
    eligible, mean, scale, rotation = model['eligible'], model['mean'], model['scale'], model['rotation']
    return np.asarray(normalized_rows[:, eligible].multiply(1 / scale) @ rotation) - (mean / scale) @ rotation


def place(normalized_rows, model):
    """No refit: frozen transform -> frozen tree projection -> LOCKED off-manifold threshold and bins from the registry (never recomputed)."""
    pcs = apply_frozen_transform(normalized_rows, model)
    p = projection(pcs, model['tree'])
    p['off_manifold'] = p['distance'] > model['threshold']
    p['state_id'] = assign_states(p, model['entries'], np.asarray(model['bins']))
    return pcs, p


def normalize_selected(cfg, features, mice, pops):
    _release_paths.expand('The S06 stage A pass (${DATA_ROOT}/helpers/s06_build.py lines 64-79) restricted to the rows of the chosen mice and populations; chunks without such rows are skipped, every other statement is the same.')
    manifest = pd.read_csv(cfg['manifest'], sep='\t', dtype={'animal_id': str}, low_memory=False, usecols=['obs_name', 'animal_id', 'cell_type_raw', 'included']).set_index('obs_name')
    pieces, rows = [], []
    with h5py.File(cfg['h5ad'], 'r') as f:
        obs = read_elem(f['obs']); var = read_elem(f['var']); m = manifest.loc[obs.index]
        keep = (m.animal_id.isin(mice) & m.cell_type_raw.isin(pops)).to_numpy()
        assert not m.included.astype(str).str.lower().eq('true').to_numpy()[keep].any(), 'a selected nucleus is an included (24-mouse) nucleus'
        positions = features.var_position.to_numpy(dtype=int)
        assert np.array_equal(var.iloc[positions].gene_ids.to_numpy(), features.stable_id.to_numpy())
        gene_lookup = np.full(len(var), -1, dtype=int); gene_lookup[positions] = np.arange(len(positions))
        g = f['layers/counts']; ptr = g['indptr'][...].astype(np.int64)
        for s, e in plan_chunks(ptr, 2000, 2_000_000):
            mask = keep[s:e]
            if not mask.any(): continue
            a, b = int(ptr[s]), int(ptr[e]); value = g['data'][a:b].astype(float); col = g['indices'][a:b]
            assert np.isfinite(value).all() and (value >= 0).all() and (value == np.floor(value)).all()
            row = np.repeat(np.arange(e - s), np.diff(ptr[s:e + 1])); totals = np.bincount(row, weights=value, minlength=e - s)
            retained = np.flatnonzero(mask); newrow = np.cumsum(mask) - 1; selected = mask[row] & (gene_lookup[col] >= 0)
            counts = sp.csr_matrix((value[selected], (newrow[row[selected]], gene_lookup[col[selected]])), shape=(len(retained), len(positions)))
            pieces.append(normalize(counts, totals[retained]))
            sub = m.iloc[s:e].iloc[retained]
            rows.append(pd.DataFrame(dict(obs_name=sub.index, animal_id=sub.animal_id.to_numpy(), cell_type_raw=sub.cell_type_raw.to_numpy(), full_count_total=totals[retained])))
    cells = pd.concat(rows, ignore_index=True); normalized = sp.vstack(pieces, format='csr')
    expected = manifest.loc[manifest.animal_id.isin(mice) & manifest.cell_type_raw.isin(pops)]
    assert len(cells) == len(expected) and normalized.shape == (len(cells), len(positions)) and set(cells.obs_name) == set(expected.index)
    return normalized, cells


def pc_means(pcs, cells, n_floor=20):
    d = pd.DataFrame(pcs, columns=[f'PC{i + 1}' for i in range(pcs.shape[1])]); d['animal_id'] = cells.animal_id.to_numpy()
    out = d.groupby('animal_id', sort=True).mean().reset_index(); out.insert(1, 'n_nuclei', d.groupby('animal_id', sort=True).size().to_numpy())
    out.insert(2, 'below_20_nuclei_flag', out.n_nuclei < n_floor)
    return out


def map_rows(pc_obs, stage_obs):
    """Stage A row of every PC row, found by obs_name through CORE/<name>/nucleus_metadata.tsv.gz; no row order is assumed."""
    lookup = pd.Series(np.arange(len(stage_obs)), index=pd.Index(stage_obs))
    assert lookup.index.is_unique and pd.Index(pc_obs).isin(lookup.index).all()
    return lookup[pc_obs].to_numpy()


def compare_parity(pcs, saved_pcs, p, saved_p, saved_off, saved_state):
    ok_pcs = bool(np.allclose(pcs, saved_pcs, rtol=1e-9, atol=1e-10))
    ok_proj = all(np.allclose(p[k], saved_p[k], rtol=1e-9, atol=1e-10, equal_nan=True) for k in ('pseudotime', 'distance', 'ambiguity_ratio')) and all(np.array_equal(p[k], saved_p[k]) for k in ('branch', 'edge', 'ambiguous'))
    ok_state = bool(np.array_equal(np.asarray(p['off_manifold'], dtype=bool), np.asarray(saved_off, dtype=bool)) and np.array_equal(np.asarray(p['state_id'], dtype=object), np.asarray(saved_state, dtype=object)))
    return ok_pcs, bool(ok_proj), ok_state


def parity_selftest(cfg, n_per_pop=1000, seed=20261007):
    """Run the SAME placement code on already-processed 24-mouse nuclei and require equality with the saved frozen outputs, so a coding error cannot hide."""
    A = cfg['A']; normalized = sp.load_npz(A / 'normalized_log10k_features.npz'); cells = pd.read_csv(A / 'cell_metadata.tsv.gz', sep='\t')
    units = pd.read_csv(A / 'unit_metadata.tsv', sep='\t', dtype={'animal_id': str}); unit = cells.unit_index.to_numpy(); rng = np.random.default_rng(seed); out = {}
    for pop in POPS:
        model = load_model(pop, cfg); name = model['name']
        nm = pd.read_csv(cfg['CORE'] / name / 'nucleus_metadata.tsv.gz', sep='\t'); stage_row = map_rows(nm.obs_name.to_numpy(), cells.obs_name.to_numpy())
        saved_pcs = np.load(cfg['CORE'] / name / 'nucleus_PCs.npy', mmap_mode='r'); saved_proj = np.load(cfg['REGISTRY'].parent / f'{name}_labelblind_projection.npz')
        saved_tab = pd.read_csv(cfg['S07PROJ'] / f'{name}_projection.tsv.gz', sep='\t', usecols=['obs_name', 'off_manifold', 'state_id'])
        assert len(nm) == len(saved_pcs) == len(saved_tab) and np.array_equal(saved_tab.obs_name.to_numpy(), nm.obs_name.to_numpy())
        j = np.sort(rng.choice(len(nm), size=min(n_per_pop, len(nm)), replace=False))
        pcs, p = place(normalized[stage_row[j]], model)
        ok_pcs, ok_proj, ok_state = compare_parity(pcs, np.asarray(saved_pcs[j]), p, {k: saved_proj[k][j] for k in saved_proj.files}, saved_tab.off_manifold.to_numpy()[j], saved_tab.state_id.to_numpy()[j])
        out[pop] = dict(n_checked=int(len(j)), pcs_equal=ok_pcs, projection_equal=ok_proj, off_manifold_and_state_equal=ok_state, max_abs_pc_diff=float(np.max(np.abs(pcs - np.asarray(saved_pcs[j])))))
    out['all_pass'] = all(v['pcs_equal'] and v['projection_equal'] and v['off_manifold_and_state_equal'] for k, v in out.items() if k in POPS)
    return out


def write_outputs(out_data, out_docs, cfg, normalized, cells, features_sha, parity, pins_before):
    out_data.mkdir(parents=True); out_docs.mkdir(parents=True)
    units = cells.groupby(['animal_id', 'cell_type_raw']).size().rename('n_nuclei').reset_index(); units['unit_index'] = np.arange(len(units)); units['age_weeks'] = 12
    units.to_csv(out_data / 'unit_metadata_12w.tsv', sep='\t', index=False)                       # label-blind: no genotype, no sex yet
    cells.to_csv(out_data / 'cell_metadata_12w.tsv.gz', sep='\t', index=False, compression='gzip'); sp.save_npz(out_data / 'normalized_log10k_features_12w.npz', normalized)
    am = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}).set_index('animal_id')[['genotype', 'sex']]      # JOINED ONLY NOW, after placement
    result = {}
    for pop in POPS:
        model = load_model(pop, cfg); sel = np.flatnonzero(cells.cell_type_raw.eq(pop).to_numpy()); sub = cells.iloc[sel].reset_index(drop=True)
        pcs, p = place(normalized[sel], model); assert np.isfinite(pcs).all() and pcs.shape[0] == len(sub)
        tag = pop.replace(' ', '_'); np.save(out_data / f'{tag}_PCs_12w.npy', pcs)
        means = pc_means(pcs, sub); means = means.join(am, on='animal_id'); means.insert(1, 'role', 'PRIMARY: PC-score mean per mouse (12 W, out-of-sample projection of the frozen discovery PCA)')
        means.to_csv(out_data / f'{tag}_PC_means_by_mouse_12w.tsv', sep='\t', index=False)
        frame = sub[['obs_name', 'animal_id']].copy(); frame['age_weeks'] = 12
        for k in ('pseudotime', 'branch', 'edge', 'distance', 'ambiguity_ratio', 'ambiguous', 'off_manifold'): frame[k] = p[k]
        frame['state_id'] = p['state_id']; frame['label'] = 'out-of-sample; not comparable to 24-mouse discovery values (locked threshold %.6g)' % model['threshold']
        frame = frame.join(am, on='animal_id'); frame.to_csv(out_data / f'{tag}_projection_12w.tsv.gz', sep='\t', index=False, compression='gzip')
        fr = frame.groupby(['animal_id', 'state_id']).size().rename('count').reset_index(); fr['total'] = fr.animal_id.map(frame.groupby('animal_id').size()); fr['fraction'] = fr['count'] / fr['total']
        fr.insert(0, 'label', 'out-of-sample'); fr.to_csv(out_data / f'{tag}_state_fractions_by_mouse_12w.tsv', sep='\t', index=False)
        # 24-mouse reference PC means from the FROZEN fitted PCs (in-sample), same layout, for comparison
        A = cfg['A']; units24 = pd.read_csv(A / 'unit_metadata.tsv', sep='\t', dtype={'animal_id': str}); nm = pd.read_csv(cfg['CORE'] / model['name'] / 'nucleus_metadata.tsv.gz', sep='\t')
        ref = pc_means(np.load(cfg['CORE'] / model['name'] / 'nucleus_PCs.npy'), units24.iloc[nm.unit_index.to_numpy()].reset_index(drop=True))
        ref = ref.join(units24.drop_duplicates('animal_id').set_index('animal_id')[['genotype', 'sex', 'age_weeks']], on='animal_id'); ref.insert(1, 'role', 'REFERENCE: 24 mice, in-sample (they fitted the discovery PCA)')
        ref.to_csv(out_data / f'{tag}_PC_means_by_mouse_24_reference.tsv', sep='\t', index=False)
        result[pop] = dict(n_nuclei=int(len(sub)), n_pcs=int(pcs.shape[1]), n_mice=int(sub.animal_id.nunique()), mice_below_20=means.loc[means.below_20_nuclei_flag, 'animal_id'].tolist(), locked_threshold=model['threshold'],
                           off_manifold_fraction=float(np.mean(p['off_manifold'])), unassigned_fraction=float(np.mean(p['state_id'] == 'UNASSIGNED')))
    return result


def reconcile(cfg, cells, mice, result, pins_before):
    manifest = pd.read_csv(cfg['manifest'], sep='\t', dtype={'animal_id': str}, low_memory=False, usecols=['obs_name', 'animal_id', 'cell_type_raw'])
    want = manifest[manifest.animal_id.isin(mice) & manifest.cell_type_raw.isin(POPS)].groupby(['animal_id', 'cell_type_raw']).size()
    got = cells.groupby(['animal_id', 'cell_type_raw']).size()
    checks = dict(per_mouse_population_counts_equal_manifest=bool(want.equals(got)), total_rows=int(len(cells)), no_duplicate_nuclei=bool(not cells.obs_name.duplicated().any()),
                  pins_unchanged_after=(verify_pins() == pins_before))
    checks['all_pass'] = all(v for k, v in checks.items() if isinstance(v, bool))
    return checks


def replay(out_data, cfg=CFG, features=None):
    """Fresh-process replay: (1) every 20th saved nucleus row re-derived straight from the h5ad by a per-row path (not the chunk path) and compared exactly; (2) PCs (bit-identical) AND the projection table (floats to 1e-12, labels and flags equal) of every 20th row recomputed from the saved matrix."""
    cells = pd.read_csv(out_data / 'cell_metadata_12w.tsv.gz', sep='\t', dtype={'animal_id': str}); normalized = sp.load_npz(out_data / 'normalized_log10k_features_12w.npz').tocsr()
    features = features if features is not None else pd.read_csv(cfg['A'] / 'feature_registry.tsv', sep='\t'); idx = np.arange(0, len(cells), 20); res = dict(n_rows_checked=int(len(idx)))
    with h5py.File(cfg['h5ad'], 'r') as f:
        obs_names = read_elem(f['obs']).index; pos = pd.Series(np.arange(len(obs_names)), index=obs_names)[cells.obs_name.to_numpy()[idx]].to_numpy()
        var = read_elem(f['var']); positions = features.var_position.to_numpy(dtype=int); lookup = np.full(len(var), -1); lookup[positions] = np.arange(len(positions)); g = f['layers/counts']; ptr = g['indptr']; same = []
        for r, i in zip(pos, idx):
            a, b = int(ptr[r]), int(ptr[r + 1]); v = g['data'][a:b].astype(float); c = g['indices'][a:b]; total = v.sum(); sel = lookup[c] >= 0
            row = normalize(sp.csr_matrix((v[sel], (np.zeros(sel.sum(), int), lookup[c][sel])), shape=(1, len(positions))), [total])
            same.append(bool((row != normalized[i]).nnz == 0 and abs(total - cells.full_count_total.iloc[i]) == 0))
    res['normalization_rows_identical'] = bool(all(same))
    pc_ok = []
    for pop in POPS:
        model = load_model(pop, cfg); sel = np.flatnonzero(cells.cell_type_raw.eq(pop).to_numpy()); take = sel[::20]
        pcs, p = place(normalized[take], model); saved = np.load(out_data / f"{pop.replace(' ', '_')}_PCs_12w.npy")[::20]; assert saved.shape == pcs.shape
        tab = pd.read_csv(out_data / f"{pop.replace(' ', '_')}_projection_12w.tsv.gz", sep='\t', dtype={'animal_id': str}).iloc[::20]; assert len(tab) == len(take)
        proj_ok = all(np.allclose(p[k], tab[k].to_numpy(dtype=float), rtol=0, atol=1e-12, equal_nan=True) for k in ('pseudotime', 'distance', 'ambiguity_ratio')) and all(np.array_equal(np.asarray(p[k]), tab[k].to_numpy()) for k in ('branch', 'edge', 'ambiguous', 'off_manifold')) and bool(np.array_equal(np.asarray(p['state_id'], dtype=object), tab['state_id'].to_numpy(dtype=object)))
        pc_ok.append(bool(np.array_equal(pcs, saved)) and bool(proj_ok))
    res['pcs_and_projection_identical'] = bool(all(pc_ok)); res['all_pass'] = res['normalization_rows_identical'] and res['pcs_and_projection_identical']
    return res


def guarded(out_docs, fn):
    """Run fn; on ANY failure after the output dirs exist, exclusive-create FAILED.json (traceback) in out_docs, then re-raise (no retry, nothing deleted)."""
    try: return fn()
    except BaseException:
        import traceback
        if Path(out_docs).exists():
            with open(Path(out_docs) / 'FAILED.json', 'x') as f: f.write(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), traceback=traceback.format_exc()), indent=2) + '\n')
        raise


def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--check', action='store_true'); g.add_argument('--smoke', action='store_true'); g.add_argument('--full', action='store_true'); g.add_argument('--replay', type=Path)
    a = ap.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'): sys.exit('refusing: the checks are asserts; run without -O')
    pins = verify_pins()
    if a.check: print(json.dumps(dict(pins_ok=True, n_pins=len(pins)))); return
    if a.replay:
        res = replay(a.replay); docs = ROOT / 'docs' / a.replay.name
        with open(docs / 'REPLAY.json', 'x') as f: f.write(json.dumps(res, indent=2) + '\n')
        print(json.dumps(res)); sys.exit(0 if res['all_pass'] else 1)
    tag = 'smoke' if a.smoke else 'full'; mice = SIX[:1] if a.smoke else SIX
    out_data = ROOT / f'results/phase4_12w_{tag}_{CFG["DATE"]}_A'; out_docs = ROOT / f'docs/phase4_12w_{tag}_{CFG["DATE"]}_A'
    assert not out_data.exists() and not out_docs.exists(), 'outputs exist; never overwrite'
    t0 = time.time(); parity = parity_selftest(CFG); assert parity['all_pass'], parity       # BEFORE any 12 W number is written
    features = pd.read_csv(CFG['A'] / 'feature_registry.tsv', sep='\t'); normalized, cells = normalize_selected(CFG, features, mice, POPS)
    def body():
        result = write_outputs(out_data, out_docs, CFG, normalized, cells, pins[_release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/feature_registry.tsv')], parity, pins)
        checks = reconcile(CFG, cells, mice, result, pins); assert checks['all_pass'], checks
        (out_docs / 'PARITY_SELFTEST.json').write_text(json.dumps(parity, indent=2) + '\n')
        (out_data / 'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(out_data.iterdir()) if p.is_file()))
        (out_docs / 'RECORD.json').write_text(json.dumps(dict(mode=tag, mice=list(mice), populations=list(POPS), seconds=time.time() - t0, result=result, reconcile=checks, pins=pins, script_sha256=sha(Path(__file__)),
                                                             note='12 W placed on frozen S07 discovery models; no refit; distance/off_manifold/state are out-of-sample; PC means are the primary output'), indent=2) + '\n')
        print(json.dumps(dict(mode=tag, seconds=round(time.time() - t0, 1), reconcile=checks['all_pass'], parity=parity['all_pass'], result=result)))
    guarded(out_docs, body)


if __name__ == '__main__':
    main()

#%% debug
