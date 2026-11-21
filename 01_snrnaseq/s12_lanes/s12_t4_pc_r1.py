
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# S12 T4 parts 1 and 2 , r1: ALONG vs ACROSS the tree for 12 discovery PCs, and the gene loadings of the same PCs. DESCRIPTIVE, READ-ONLY, no interpretation; nothing counts toward Step 11 (C2).
# Populations and PCs (pre-specified): Purkinje cell PC 5, 6, 11, 30; Bergmann glia PC 2, 4, 5, 6, 8, 9, 20, 29.
# PART 1. The label-blind registry arrays (projected_point, distance), the discovery PC scores and the discovery projection file are proved row-aligned (||PC - projected_point|| per nucleus equals the stored distance to 1e-9; the npz distance equals the projection-file distance).
#   Per nucleus: on-tree part = projected_point[:, k]; off-tree part = (PC - projected_point)[:, k]. Per mouse: the mean of each part and of the full PC, checked against the T3 table. KI vs WT within age with the T3 exact test (separation, exact permutation p).
# PART 2. For each PC the 25 most positive and 25 most negative loadings (pca_model.npz rotation), mapped to genes ONLY through the rows of ${DATA_ROOT}/results/s06_stageA_v01/feature_registry.tsv after PROVING that feature_indices index those rows (see PROOF below),
#   each gene flagged HIT (with ages) / TESTED_NOT_HIT / NOT_IN_RELEASE from the S11 release 40593a55 (decision stripped and upper-cased == HIT and failed_fit_p1 false; any age; no family filter), the ambient_flagged_hit label carried, and module membership from MODULE_IDS.json (stable_id keyed). Lists only; no enrichment test.
# PROOF (STOP on any failure): (a) sha256 of S06 SHA256SUMS.txt equals source_transform_registry_sha256 in the S07 core RECORD.json, and that list holds the pinned hashes of feature_registry.tsv, transform_parameters.npz, transform_keys.json and normalized_log10k_features.npz;
#   (b) the registry has exactly 1414 rows with var_position strictly increasing; (c) the stored shape of normalized_log10k_features.npz is (280295, 1414); (d) for the (discovery, population) key of transform_keys.json, feature_indices == flatnonzero(transform_parameters['eligible'][key index])
#   and rotation has shape (len(feature_indices), 30); (e) every feature index is below 1414.
# usage: python -I -B s12_t4_pc_r1.py --real OUT   |   python -I -B s12_t4_pc_r1.py --fixture FIXDIR OUT      (OUT must exist and be empty)
import argparse, hashlib, itertools, json, os, re, sys, time, zipfile
import numpy as np
import pandas as pd

ROOT = _release_paths.expand('${REPO_ROOT}/')
POPS = [('Purkinje cell', 'Purkinje_cell', 'Purkinje_cell', [5, 6, 11, 30]), ('Bergmann glia', 'Bergmann_glia', 'Bergmann_glia', [2, 4, 5, 6, 8, 9, 20, 29])]       # (S07 / release registry name, file stem, T3 name, PCs)
AGES = [5, 18, 24, 30]
NPC, NFEAT, NCELLS, TOPN, FLOOR = 30, 1414, 280295, 25, 20
S6 = _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/')
PIN_MAN = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv')
PIN_T3 = _release_paths.expand('${REPO_ROOT}/results/s12_lanes/mouse_table.tsv')
PIN_REL = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz')
PIN_MOD = _release_paths.expand('${DATA_ROOT}/s08_module_ids/MODULE_IDS.json')
PIN_REC = _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/RECORD.json')
REAL_PINS = {
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz'): "40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298",
    _release_paths.expand('${DATA_ROOT}/s08_module_ids/MODULE_IDS.json'): "458713e72b81097063d8b2d51bd784b11b125acee614e65388a815dee15b8866",
    _release_paths.expand('${REPO_ROOT}/results/s12_lanes/mouse_table.tsv'): "94cdfb6f7213dfee8d9ad01f30c0b6f99404335007a234cb2c4dc7b6d3938b00",
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): "7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/SHA256SUMS.txt'): "79cfae1a1bbf9f90c407425636600e31bbe1d9913f61698c3473b738212bcff0",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/feature_registry.tsv'): "d49b138d3a806db64fb37b6994a30034569cd530bfe3b6a82becf10f367f8c24",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/normalized_log10k_features.npz'): "082f1639be27db8791529bcaae7710b9cca4f8bf85b3f563ddd31de9455d9507",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/transform_keys.json'): "a9f42ea00582c9e9174c732cb3c9ee9090507bbcad53cd442d23c3d5f28889b8",
    _release_paths.expand('${DATA_ROOT}/results/s06_stageA_v01/transform_parameters.npz'): "0e4b000d17a09e33b9277f8e19e00eabf9e8dab319cbd71e179ad4b42cdb3b3b",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/RECORD.json'): "594459e0e5a1b3a30f2f7b8e284968036aae3fafbaf02647c06871d9c5ec2a2a",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_PCs.npy'): "ed535ac9cdb36c51498155641172de92d9b6daa65e611ffa81a705b1a573ea13",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): "23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/pca_model.npz'): "34c1e8e2074b3c8247b09fac35d64acd6070e9f745b6604efdac3175f9642596",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'): "f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): "0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/pca_model.npz'): "eaf1468f99ed7935acf94a1e58b1a45f1fe1333ae6cc541c6afc550e4d128284",
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Bergmann_glia_projection.tsv.gz'): "29b779ee72e327e066e5e9efa8ece185560498b1dbca650cd71aa9019342ae10",
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Purkinje_cell_projection.tsv.gz'): "8605f82700db2e783032ce66dd744be479af987186595a0d6d09d4a860c6b458",
    _release_paths.expand('${DATA_ROOT}/results/s07_state_registry_v01/discovery__Bergmann_glia_labelblind_projection.npz'): "71534ccdbdced588783e1140b89544a733df9e35d6264087bd63284ab712ecf4",
    _release_paths.expand('${DATA_ROOT}/results/s07_state_registry_v01/discovery__Purkinje_cell_labelblind_projection.npz'): "219ea6843996fc4a2f0f863c7a653efe7e51cbf6ce0bc30e6798050e749bb00f",
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


def r(rel_, kind, stem=None):
    return {'npz': _release_paths.expand('${DATA_ROOT}/results/s07_state_registry_v01/discovery__%s_labelblind_projection.npz'), 'pcs': _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__%s/nucleus_PCs.npy'), 'meta': _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__%s/nucleus_metadata.tsv.gz'),
            'disc': _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__%s_projection.tsv.gz'), 'pca': _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__%s/pca_model.npz')}[kind] % stem


def npz_shapes(p):
    out = {}
    with zipfile.ZipFile(p) as z:
        for m in z.namelist():
            if m.endswith('.npy'):
                with z.open(m) as f:
                    v = np.lib.format.read_magic(f); sh, fo, dt = (np.lib.format.read_array_header_1_0(f) if v == (1, 0) else np.lib.format.read_array_header_2_0(f)); out[m[:-4]] = tuple(sh)
    return out


def geno(s):
    k = str(s).strip().upper()
    if k in ('WT', 'KI'): return k
    stop('unrecognised genotype label: ' + repr(s))


def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True); g.add_argument('--real', action='store_true'); g.add_argument('--fixture'); ap.add_argument('out')
    a = ap.parse_args(); t0 = time.time()
    if a.real: base, pins = '', REAL_PINS
    else: base = a.fixture.rstrip('/') + '/'; pins = json.load(open(base + 'FIXTURE_PINS.json'))
    if not os.path.isdir(a.out) or os.listdir(a.out): stop('the output directory must exist and be empty')
    out = a.out.rstrip('/') + '/'
    for p, h in pins.items():
        if not os.path.isfile(base + p): stop('missing input ' + p)
        if sha256(base + p) != h: stop('sha256 differs from the pin: ' + p)
    want = [r(None, k, s) for _, s, _, _ in POPS for k in ('npz', 'pcs', 'meta', 'disc', 'pca')] + [S6 + x for x in ('feature_registry.tsv', 'transform_parameters.npz', 'transform_keys.json', 'normalized_log10k_features.npz', 'SHA256SUMS.txt')] + [PIN_MAN, PIN_T3, PIN_REL, PIN_MOD, PIN_REC]
    if sorted(want) != sorted(pins): stop('the pin list is not exactly the 20 inputs of this script')
    ok('all %d inputs match their pinned sha256' % len(pins), len(pins))
    # ---------------- manifest and T3 table
    man = pd.read_csv(base + PIN_MAN, sep='\t', dtype=str, keep_default_na=False)
    man = man[man['included'].str.strip().str.lower().isin(['true', '1', 'yes', 'y', 't'])].copy(); man['animal_id'] = man['animal_id'].str.strip(); man['genotype'] = man['genotype'].map(geno); man['age_weeks'] = man['age_weeks'].astype(float).astype(int)
    if len(man) != 24 or man['animal_id'].duplicated().any() or not (man.groupby(['age_weeks', 'genotype']).size() == 3).all(): stop('the pinned manifest is not 24 unique mice, 3 WT + 3 KI per age')
    man = man.set_index('animal_id'); ok('pinned manifest: 24 included mice, 3 WT + 3 KI per age', 24)
    t3 = pd.read_csv(base + PIN_T3, sep='\t', dtype={'animal_id': str})
    # ---------------- PART 1
    mean_rows = []; test_rows = []; part1 = {}
    for name, stem, t3name, pcs_sel in POPS:
        z = np.load(base + r(None, 'npz', stem)); D = pd.read_csv(base + r(None, 'disc', stem), sep='\t', dtype={'animal_id': str}); meta = pd.read_csv(base + r(None, 'meta', stem), sep='\t'); pc = np.load(base + r(None, 'pcs', stem))
        for need in ('projected_point', 'distance'):
            if need not in z.files: stop('%s: the registry npz lacks %s' % (name, need))
        pp, dz = z['projected_point'], z['distance']; n = len(D)
        if pc.shape != (n, NPC) or pp.shape != (n, NPC) or dz.shape != (n,) or len(meta) != n: stop(name + ': row counts or shapes differ between the registry npz, nucleus_PCs, nucleus_metadata and the discovery projection')
        if meta['obs_name'].tolist() != D['obs_name'].tolist() or D['obs_name'].duplicated().any(): stop(name + ': obs_name order differs between nucleus_metadata and the discovery projection')
        if not np.allclose(dz, D['distance'].to_numpy(float), rtol=0, atol=1e-12): stop(name + ': the registry npz distance differs from the discovery projection distance')
        if not np.allclose(np.linalg.norm(pc - pp, axis=1), dz, rtol=0, atol=1e-9): stop(name + ': ||PC - projected_point|| differs from the stored distance (rows are not aligned)')
        ok(name + ': nucleus_PCs and projected_point are row-aligned (||PC - projected_point|| equals the stored distance to 1e-9; npz distance equals the projection-file distance)', n)
        first = D.groupby('animal_id')[['age_weeks', 'genotype']].first()
        if not (D.groupby('animal_id')[['age_weeks', 'genotype']].nunique() == 1).all().all() or not set(first.index) <= set(man.index): stop(name + ': a mouse has more than one age or genotype, or is not an included manifest mouse')
        for aid, rr in first.iterrows():
            if geno(rr['genotype']) != man.loc[aid, 'genotype'] or int(rr['age_weeks']) != int(man.loc[aid, 'age_weeks']): stop('%s: mouse %s disagrees with the pinned manifest' % (name, aid))
        ok(name + ': genotype and age of every mouse equal the pinned manifest', len(first))
        T = t3[t3['population'] == t3name].set_index('animal_id')
        if set(T.index) != set(first.index): stop(name + ': the T3 table and the projection file hold different mice')
        for aid, rr in first.iterrows():
            sel = (D['animal_id'] == aid).to_numpy(); nn = int(sel.sum())
            for k in pcs_sel:
                full = float(pc[sel, k - 1].mean()); on = float(pp[sel, k - 1].mean()); off = float((pc[sel, k - 1] - pp[sel, k - 1]).mean())
                if abs(full - (on + off)) > 1e-9: stop('%s %s PC%d: full mean differs from on-tree plus off-tree mean' % (name, aid, k))
                t3v = float(T.loc[aid, 'discovery_PC%d_mean' % k])
                if abs(t3v - full) > 1e-8 * max(1.0, abs(t3v)): stop('%s %s PC%d: full PC mean differs from the T3 table value' % (name, aid, k))
                mean_rows.append(dict(population=name, animal_id=aid, genotype=man.loc[aid, 'genotype'], age_weeks=int(man.loc[aid, 'age_weeks']), n_nuclei=nn, included=nn >= FLOOR, PC=k, full_mean=full, on_tree_mean=on, off_tree_mean=off))
        ok(name + ': per-mouse full PC means equal the T3 table; full = on-tree + off-tree', len(first) * len(pcs_sel))
    MM = pd.DataFrame(mean_rows)
    for name, stem, t3name, pcs_sel in POPS:
        for k in pcs_sel:
            for age in AGES:
                sub = MM[(MM['population'] == name) & (MM['PC'] == k) & (MM['age_weeks'] == age) & MM['included']]
                for part in ('full_mean', 'on_tree_mean', 'off_tree_mean'):
                    w = sub[sub['genotype'] == 'WT']; ki = sub[sub['genotype'] == 'KI']; wv = w[part].to_numpy(float); kv = ki[part].to_numpy(float)
                    rr = dict(population=name, PC=k, part=part.replace('_mean', ''), age_weeks=age, n_WT=len(w), n_KI=len(ki))
                    if len(w) >= 1 and len(ki) >= 1:
                        rr.update(WT_min=wv.min(), WT_max=wv.max(), WT_mean=wv.mean(), KI_min=kv.min(), KI_max=kv.max(), KI_mean=kv.mean(), diff_KI_minus_WT=kv.mean() - wv.mean(), separation=separation(wv, kv),
                                  WT_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(w['animal_id'], wv)), KI_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(ki['animal_id'], kv)))
                        pr = exact_perm(wv, kv)
                        if pr: rr.update(pr)
                    test_rows.append(rr)
    TT = pd.DataFrame(test_rows)
    summ = TT.pivot_table(index=['population', 'PC', 'age_weeks'], columns='part', values='separation', aggfunc='first').reset_index(); summ.columns = ['sep_' + c if c in ('full', 'on_tree', 'off_tree') else c for c in summ.columns]
    # ---------------- PART 2: proof
    sums = open(base + S6 + 'SHA256SUMS.txt').read().splitlines(); listed = {l.split(maxsplit=1)[1].strip(): l.split()[0] for l in sums if l.strip()}
    rec = json.load(open(base + PIN_REC))
    if sha256(base + S6 + 'SHA256SUMS.txt') != rec.get('source_transform_registry_sha256'): stop('S06 SHA256SUMS.txt does not equal source_transform_registry_sha256 of the S07 core RECORD.json')
    for f in ('feature_registry.tsv', 'transform_parameters.npz', 'transform_keys.json', 'normalized_log10k_features.npz'):
        if listed.get(f) != pins[S6 + f]: stop('S06 SHA256SUMS.txt does not list the pinned hash of ' + f)
    ok('proof (a): the S07 core RECORD pins the S06 stage A checksum list, which lists the pinned hashes of the registry, transform parameters, transform keys and normalized features', 5)
    reg = pd.read_csv(base + S6 + 'feature_registry.tsv', sep='\t', dtype={'var_position': int})
    if len(reg) != NFEAT or not (np.diff(reg['var_position'].to_numpy()) > 0).all(): stop('proof (b): feature_registry.tsv does not have exactly %d rows with strictly increasing var_position' % NFEAT)
    ok('proof (b): feature_registry.tsv has %d rows, var_position strictly increasing' % NFEAT, NFEAT)
    with zipfile.ZipFile(base + S6 + 'normalized_log10k_features.npz') as z:
        with z.open('shape.npy') as f: shp = tuple(int(x) for x in np.lib.format.read_array(f))
    if shp != (NCELLS, NFEAT): stop('proof (c): normalized_log10k_features.npz shape %s is not (%d, %d)' % (shp, NCELLS, NFEAT))
    ok('proof (c): normalized_log10k_features.npz stored shape is (%d, %d)' % (NCELLS, NFEAT), 1)
    keys = json.load(open(base + S6 + 'transform_keys.json')); tp = np.load(base + S6 + 'transform_parameters.npz'); elig = tp['eligible']
    if elig.shape != (len(keys), NFEAT): stop('proof (d): eligible has shape %s for %d keys' % (elig.shape, len(keys)))
    rel = pd.read_csv(base + PIN_REL, sep='\t', dtype=str, keep_default_na=False, usecols=['stable_id', 'population_registry_name', 'age_weeks', 'family', 'decision', 'failed_fit_p1', 'ambient_flagged_hit', 'dataset_symbol'])
    mod = json.load(open(base + PIN_MOD))['modules']
    for mname, e in mod.items():
        if not isinstance(e.get('stable_ids'), list) or not all(isinstance(s, str) and s.startswith('ENSMUSG') for s in e['stable_ids']): stop('MODULE_IDS.json: module %s is not a list of stable_id strings' % mname)
    ok('MODULE_IDS.json: every module holds a list of ENSMUSG stable_id strings', len(mod))
    rel['dec'] = rel['decision'].str.strip().str.upper(); rel['failed'] = rel['failed_fit_p1'].str.strip().str.lower()
    if not set(rel['failed']) <= {'false', 'true'}: stop('unexpected failed_fit_p1 labels: ' + ', '.join(sorted(set(rel['failed']))))
    rel['is_hit'] = (rel['dec'] == 'HIT') & (rel['failed'] == 'false'); rel['age'] = rel['age_weeks'].astype(float).astype(int)
    amb_labels = sorted(rel['ambient_flagged_hit'].unique())
    modsets = {m: set(e['stable_ids']) for m, e in mod.items()}; load_rows = []; align_rows = []
    for name, stem, t3name, pcs_sel in POPS:
        pm = np.load(base + r(None, 'pca', stem)); fi, rot = pm['feature_indices'], pm['rotation']
        ks = [i for i, kk in enumerate(keys) if kk['mode'] == 'discovery' and kk['population'] == name]
        if len(ks) != 1: stop('proof (d): %s has %d discovery keys in transform_keys.json' % (name, len(ks)))
        if not np.array_equal(fi, np.flatnonzero(elig[ks[0]])) or rot.shape != (len(fi), NPC) or fi.min() < 0 or fi.max() >= NFEAT: stop('proof (d/e): %s feature_indices differ from flatnonzero(eligible[key]) or rotation has the wrong shape' % name)
        ok('proof (d, e): %s feature_indices equal flatnonzero(eligible[%d]) (%d features) and rotation is (%d, %d)' % (name, ks[0], len(fi), len(fi), NPC), len(fi))
        align_rows.append(dict(population=name, transform_key_index=ks[0], n_features=len(fi), rotation_shape='%dx%d' % rot.shape, feature_index_min=int(fi.min()), feature_index_max=int(fi.max())))
        R_ = rel[rel['population_registry_name'] == name]
        if R_.empty: stop(name + ': no rows of this population in the S11 release')
        for k in pcs_sel:
            col = rot[:, k - 1]; order = np.argsort(col, kind='stable'); neg = order[:TOPN]; pos = order[::-1][:TOPN]
            for direction, idxs in (('positive', pos), ('negative', neg)):
                for rank, j in enumerate(idxs, 1):
                    fx = int(fi[j]); row = reg.iloc[fx]; sid = row['stable_id']; g_ = R_[R_['stable_id'] == sid]
                    if g_.empty: state, ages, amb = 'NOT_IN_RELEASE', '', ''
                    else:
                        h_ = g_[g_['is_hit']]; state = 'HIT' if len(h_) else 'TESTED_NOT_HIT'; ages = ';'.join(str(x) for x in sorted(h_['age'])); amb = ';'.join(sorted(set(h_['ambient_flagged_hit']))) if len(h_) else ''
                    mods = ';'.join(m for m in modsets if sid in modsets[m])
                    load_rows.append(dict(population=name, PC=k, direction=direction, rank=rank, loading=float(col[j]), feature_index=fx, stable_id=sid, dataset_symbol=row['dataset_symbol'], canonical_symbols=row['canonical_symbols'], core_symbols=row['core_symbols'], s08_state=state, s08_hit_ages=ages, ambient_flagged_hit_labels_on_hit_rows=amb, modules=mods))
    LD = pd.DataFrame(load_rows)
    def wr(df, nm): df.to_csv(out + nm, sep='\t', index=False, float_format='%.10g', mode='x')
    wr(MM, 'along_across_mouse_means.tsv'); wr(TT, 'along_across_kivswt.tsv'); wr(summ, 'along_across_summary.tsv'); wr(LD, 'top_loadings.tsv'); wr(pd.DataFrame(align_rows), 'alignment_proof.tsv'); wr(pd.DataFrame(CHECKS, columns=['check_passed', 'items']), 'checks.tsv')
    res = dict(mode='real' if a.real else 'fixture', checks=len(CHECKS), release_labels=dict(population_registry_name=sorted(rel['population_registry_name'].unique()), family=sorted(rel['family'].unique()), decision=sorted(rel['decision'].unique()), failed_fit_p1=sorted(rel['failed_fit_p1'].unique()), ambient_flagged_hit=amb_labels),
               hit_rule='decision stripped and upper-cased == HIT and failed_fit_p1 false; any age; no family filter', modules_in_file=sorted(mod), tests=len(TT), loading_rows=len(LD), s08_states={s: int((LD['s08_state'] == s).sum()) for s in ('HIT', 'TESTED_NOT_HIT', 'NOT_IN_RELEASE')},
               note='lists and counts only; no enrichment test; descriptive, not counted (S11 = C2)', python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, script_sha256=sha256(os.path.abspath(__file__)), seconds=round(time.time() - t0, 2))
    json.dump(res, open(out + 'RESULT.json', 'x'), indent=1, sort_keys=True); print(out)


if __name__ == '__main__': main()
#%% debug
