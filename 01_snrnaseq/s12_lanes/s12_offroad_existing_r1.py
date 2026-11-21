
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# S12.2 existing-output test , r1. DESCRIPTIVE; nothing here counts toward Step 11 (C2).
# Purkinje first, then Bergmann glia. Per mouse (24 mice, no 12 W): median HELD-OUT distance to the principal tree (fold files, rows with training == False), the same divided by the fold's stored distance_threshold (a proposal, because each fold has its own PC space and tree),
# held-out off_manifold_fraction, discovery-file distance and off_manifold_fraction (secondary: the discovery model was fitted on every nucleus), nuclei counts, and the per-mouse mean of each of the 30 discovery PC scores.
# KI vs WT within age (never pooled across ages): per-arm ranges, complete separation (yes/no and direction), exact permutation p over all labelings of the mice (20 for 3 v 3). No projection, fit, residual or nucleus-level test is computed.
# Inputs are only the pinned S07 files (full sha256, checked before any file is parsed). Every STOP below is final: no retry. Outputs are written only into a NEW directory (exclusive create).
# usage: python -I -B s12_offroad_existing_r1.py --real OUT   |   python -I -B s12_offroad_existing_r1.py --fixture FIXDIR OUT      (FIXDIR holds FIXTURE_PINS.json and the same relative file layout)
import argparse, hashlib, itertools, json, math, os, sys, time
import numpy as np
import pandas as pd

ROOT = _release_paths.expand('${REPO_ROOT}/')
POPS = ['Purkinje_cell', 'Bergmann_glia']
AGES = [5, 18, 24, 30]
FLOOR = 20                       # PROPOSAL: minimum held-out nuclei per mouse and population
NPC = 30
REAL_PINS = {
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): "7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_PCs.npy'): "ed535ac9cdb36c51498155641172de92d9b6daa65e611ffa81a705b1a573ea13",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): "23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'): "f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b",
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): "0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721",
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Bergmann_glia_projection.tsv.gz'): "29b779ee72e327e066e5e9efa8ece185560498b1dbca650cd71aa9019342ae10",
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__Purkinje_cell_projection.tsv.gz'): "8605f82700db2e783032ce66dd744be479af987186595a0d6d09d4a860c6b458",
    _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/projection_mouse_diagnostics.tsv'): "59801b69c6a0fe9d46f036d8f74d45dba20c0c60d1535830f355ba5b38c686ac",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold1__Bergmann_glia_projection.tsv.gz'): "a66b0ddb1e911d4b84658cecf9aa61a95520e29a066ac025dc5150aede6119c6",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold1__Purkinje_cell_projection.tsv.gz'): "171e36bb9e0dae04ba2d3c42007d8d414cf753b4e02a6e58455280e49542a9dc",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold2__Bergmann_glia_projection.tsv.gz'): "ba844daebff702653cb6ef8118d94d5dc8b81585ff936f5629446675ea9e25ac",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold2__Purkinje_cell_projection.tsv.gz'): "d5e92bc563488c57759de94c15307d34fb01a04d50d948ddef4f69f76cb40f96",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold3__Bergmann_glia_projection.tsv.gz'): "a1cc63913a6939f2a9def942d62d33d0949746472a6d93f7063ca6413324cd21",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold3__Purkinje_cell_projection.tsv.gz'): "28df4f33d26f3646af99deef046fb05f48b1ada13a64c14c9a601b1b65399289",
    _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/projection_mouse_diagnostics.tsv'): "c5d5bffc54f019e9c9055ad592fe31a1590731f4af51b76ba569d01f84a934d3",
}
CHECKS = []


def stop(msg): raise SystemExit('STOP: ' + msg)


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()


def rel(kind, pop=None, k=None):
    S = {'fold': _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/fold%d__%s_projection.tsv.gz'), 'disc': _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/discovery__%s_projection.tsv.gz'),
         'pcs': _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__%s/nucleus_PCs.npy'), 'meta': _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__%s/nucleus_metadata.tsv.gz')}
    if kind == 'fold': return S['fold'] % (k, pop)
    if kind in S: return S[kind] % pop
    return {'diag_disc': _release_paths.expand('${DATA_ROOT}/results/s07_discovery_projection_v02/projection_mouse_diagnostics.tsv'), 'diag_fold': _release_paths.expand('${DATA_ROOT}/results/s07_fold_projection_v01/projection_mouse_diagnostics.tsv'),
            'manifest': _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv')}[kind]


def exact_perm(wt, ki):
    """All labelings of the mice into arms of the observed sizes; statistic = mean(KI) - mean(WT); ties within 1e-12 count as equal."""
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


def read_tsv(p, **kw): return pd.read_csv(p, sep='\t', **kw)


def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True); g.add_argument('--real', action='store_true'); g.add_argument('--fixture'); ap.add_argument('out')
    a = ap.parse_args(); t0 = time.time()
    if a.real: base, pins = '', REAL_PINS
    else:
        base = a.fixture.rstrip('/') + '/'; pins = json.load(open(base + 'FIXTURE_PINS.json'))
    if not os.path.isdir(a.out) or os.listdir(a.out): stop('the output directory must exist and be empty')
    out = a.out.rstrip('/') + '/'
    # ---- pins: every input file is hashed and compared BEFORE any file is parsed
    for p, h in pins.items():
        if not os.path.isfile(base + p): stop('missing input ' + p)
        if sha256(base + p) != h: stop('sha256 differs from the pin: ' + p)
    wanted = [rel('fold', pop, k) for pop in POPS for k in (1, 2, 3)] + [rel(x, pop) for pop in POPS for x in ('disc', 'pcs', 'meta')] + [rel('diag_disc'), rel('diag_fold'), rel('manifest')]
    if sorted(wanted) != sorted(pins): stop('the pin list is not exactly the 15 inputs of this script')
    CHECKS.append(('all %d inputs match their pinned sha256' % len(pins), len(pins)))
    # ---- manifest (pinned animal manifest)
    m = read_tsv(base + rel('manifest'), dtype=str, keep_default_na=False)
    for c in ['animal_id', 'genotype', 'age_weeks', 'sex', 'included']:
        if c not in m.columns: stop('manifest lacks column ' + c)
    m = m[m['included'].str.strip().str.lower().isin(['true', '1', 'yes', 'y', 't'])].copy()
    m['animal_id'] = m['animal_id'].str.strip(); m['genotype'] = m['genotype'].str.strip().str.upper(); m['sex'] = m['sex'].str.strip().str.upper().str[0]; m['age_weeks'] = m['age_weeks'].astype(float).astype(int)
    CHECKS.append(('manifest: 24 included mice', 1)) if len(m) == 24 else stop('manifest has %d included mice, expected 24' % len(m))
    if set(m['genotype']) != {'WT', 'KI'} or set(m['age_weeks']) != set(AGES): stop('manifest genotypes or ages are not WT/KI and 5/18/24/30')
    if not (m.groupby(['age_weeks', 'genotype']).size() == 3).all() or m['animal_id'].duplicated().any(): stop('manifest is not 3 WT + 3 KI per age with unique animal ids')
    CHECKS.append(('manifest: 3 WT + 3 KI at each of 5, 18, 24, 30 W; unique ids', 24))
    man = m.set_index('animal_id')
    # ---- diagnostics tables
    dd = read_tsv(base + rel('diag_disc'), dtype={'animal_id': str}); df = read_tsv(base + rel('diag_fold'), dtype={'animal_id': str})
    for t in (dd, df): t['pop'] = t['population'].str.replace(' ', '_')
    fold_ofm = df.set_index(['mode', 'pop', 'animal_id'])['off_manifold_fraction']; disc_ofm = dd.set_index(['mode', 'pop', 'animal_id'])['off_manifold_fraction']      # stored per-mouse values
    mouse_rows = []; test_rows = []; pcs_store = {}
    for pop in POPS:
        D = read_tsv(base + rel('disc', pop), dtype={'animal_id': str}); F = {k: read_tsv(base + rel('fold', pop, k), dtype={'animal_id': str}) for k in (1, 2, 3)}
        meta = read_tsv(base + rel('meta', pop)); pcs = np.load(base + rel('pcs', pop), mmap_mode='r')
        n = len(D)
        if D['obs_name'].duplicated().any(): stop(pop + ': duplicate obs_name in the discovery file')
        if len(meta) != n or pcs.shape != (n, NPC): stop(pop + ': discovery PCs or metadata do not have the rows of the discovery projection')
        if meta['obs_name'].tolist() != D['obs_name'].tolist(): stop(pop + ': obs_name order differs between nucleus_metadata and the discovery projection')
        CHECKS.append((pop + ': discovery projection, nucleus_metadata and nucleus_PCs rows are aligned (obs_name, shape n x 30)', n))
        for k in (1, 2, 3):
            if F[k]['obs_name'].tolist() != D['obs_name'].tolist(): stop('%s fold%d: obs_name set or order differs from the discovery file' % (pop, k))
            for c in ['animal_id', 'age_weeks', 'genotype', 'sex', 'fold']:
                if not (F[k][c].astype(str).to_numpy() == D[c].astype(str).to_numpy()).all(): stop('%s fold%d: column %s differs from the discovery file' % (pop, k, c))
        CHECKS.append((pop + ': the three fold files carry the same nuclei, order and mouse columns as the discovery file', 3 * n))
        # per-mouse attributes: unique within a mouse, equal to the manifest
        att = D.groupby('animal_id')[['age_weeks', 'genotype', 'sex', 'fold']].nunique()
        if not (att == 1).all().all(): stop(pop + ': a mouse has more than one age, genotype, sex or fold')
        first = D.groupby('animal_id')[['age_weeks', 'genotype', 'sex', 'fold']].first()
        if not set(first.index) <= set(man.index): stop(pop + ': a mouse in the projection files is not an included mouse of the pinned manifest')
        for aid, r in first.iterrows():
            if str(r['genotype']).strip().upper() != man.loc[aid, 'genotype'] or int(r['age_weeks']) != int(man.loc[aid, 'age_weeks']) or str(r['sex']).strip().upper()[0] != man.loc[aid, 'sex']: stop('%s: mouse %s disagrees with the pinned manifest in genotype, age or sex' % (pop, aid))
        CHECKS.append((pop + ': genotype, age and sex of every mouse equal the pinned manifest', len(first)))
        if not first['fold'].isin([1, 2, 3]).all(): stop(pop + ': a mouse has a fold outside 1, 2, 3')
        # held out exactly once
        for k in (1, 2, 3):
            if not (F[k]['training'].to_numpy() == (F[k]['fold'].to_numpy() != k)).all(): stop('%s fold%d: training flag differs from (fold != %d)' % (pop, k, k))
        held = sum((~F[k]['training'].to_numpy()).astype(int) for k in (1, 2, 3))
        if not (held == 1).all(): stop(pop + ': a nucleus is not held out exactly once across the three fold files')
        CHECKS.append((pop + ': every nucleus is held out in exactly one fold file', n))
        if not D['training'].all(): stop(pop + ': the discovery file has a nucleus with training == False')
        for t, name, f in [(D, 'discovery', None)] + [(F[k], 'fold%d' % k, k) for k in (1, 2, 3)]:
            if not np.isfinite(t['distance'].to_numpy()).all(): stop('%s %s: non-finite distance' % (pop, name))
        # stored thresholds and diagnostics reproduce the stored flags
        thr = {}
        for name, t, dg in [('discovery', D, dd[(dd['mode'] == 'discovery') & (dd['pop'] == pop)])] + [('fold%d' % k, F[k], df[(df['mode'] == 'fold%d' % k) & (df['pop'] == pop)]) for k in (1, 2, 3)]:
            if dg.empty or dg['distance_threshold'].nunique() != 1: stop('%s %s: diagnostics rows missing or more than one threshold' % (pop, name))
            th = float(dg['distance_threshold'].iloc[0]); thr[name] = th
            if not (t['off_manifold'].to_numpy() == (t['distance'].to_numpy() > th)).all(): stop('%s %s: stored off_manifold differs from distance > the stored threshold' % (pop, name))
            g = t.groupby('animal_id'); frac = g['off_manifold'].mean(); cnt = g.size(); trn = g['training'].first(); dgi = dg.set_index('animal_id')
            if set(dgi.index) != set(frac.index): stop('%s %s: diagnostics mice differ from the projection file mice' % (pop, name))
            if not ((dgi.loc[frac.index, 'n_nuclei'].to_numpy() == cnt.to_numpy()).all() and (dgi.loc[frac.index, 'training'].to_numpy() == trn.to_numpy()).all() and np.allclose(dgi.loc[frac.index, 'off_manifold_fraction'].to_numpy(), frac.to_numpy(), rtol=0, atol=1e-12)): stop('%s %s: diagnostics n_nuclei, training or off_manifold_fraction differ from the projection file' % (pop, name))
        CHECKS.append((pop + ': stored off_manifold equals distance > stored threshold, and the diagnostics rows reproduce n_nuclei, training and off_manifold_fraction, in the discovery and all three fold files', 4 * n))
        # per-mouse measures
        pc = np.asarray(pcs)
        for aid, r in first.iterrows():
            kfold = int(r['fold']); sel_d = (D['animal_id'] == aid).to_numpy(); Fk = F[kfold]; sel_h = ((Fk['animal_id'] == aid) & ~Fk['training']).to_numpy()
            if sel_h.sum() != sel_d.sum(): stop('%s %s: held-out rows differ from the mouse nuclei' % (pop, aid))
            n_h = int(sel_h.sum()); incl = n_h >= FLOOR
            med = float(np.median(Fk['distance'].to_numpy()[sel_h])); frac_h = float(fold_ofm[('fold%d' % kfold, pop, aid)])        # diagnostics row: mode fold<own fold>, training False
            row = dict(population=pop, animal_id=aid, genotype=man.loc[aid, 'genotype'], age_weeks=int(man.loc[aid, 'age_weeks']), sex=man.loc[aid, 'sex'], fold=kfold, n_nuclei=int(sel_d.sum()), n_heldout=n_h, included=incl,
                       heldout_median_distance=med, fold_distance_threshold=thr['fold%d' % kfold], heldout_median_distance_over_threshold=med / thr['fold%d' % kfold], heldout_off_manifold_fraction=frac_h,
                       discovery_median_distance=float(np.median(D['distance'].to_numpy()[sel_d])), discovery_off_manifold_fraction=float(disc_ofm[('discovery', pop, aid)]), diag_mode_used='fold%d (heldout) and discovery' % kfold)
            for j in range(NPC): row['discovery_PC%d_mean' % (j + 1)] = float(pc[sel_d, j].mean())
            mouse_rows.append(row)
    MT = pd.DataFrame(mouse_rows)
    # ---- KI vs WT within age
    measures = [('heldout_median_distance', 'primary'), ('heldout_median_distance_over_threshold', 'primary_normalised_PROPOSAL'), ('heldout_off_manifold_fraction', 'primary_second_measure'),
                ('discovery_median_distance', 'secondary_not_held_out'), ('discovery_off_manifold_fraction', 'secondary_not_held_out'), ('n_nuclei', 'info_capture'), ('n_heldout', 'info_capture')] + [('discovery_PC%d_mean' % (j + 1), 'PC_check') for j in range(NPC)]
    for pop in POPS:
        for age in AGES:
            sub = MT[(MT['population'] == pop) & (MT['age_weeks'] == age) & MT['included']]
            for mname, role in measures:
                w = sub[sub['genotype'] == 'WT']; ki = sub[sub['genotype'] == 'KI']; wv = w[mname].to_numpy(float); kv = ki[mname].to_numpy(float)
                r = dict(population=pop, age_weeks=age, measure=mname, role=role, family={'primary': 'primary_24', 'primary_normalised_PROPOSAL': 'primary_24', 'primary_second_measure': 'primary_24', 'PC_check': 'PC_check_240'}.get(role, 'other'), n_WT=len(w), n_KI=len(ki))
                if len(w) >= 1 and len(ki) >= 1:
                    r.update(WT_min=wv.min(), WT_max=wv.max(), WT_mean=wv.mean(), KI_min=kv.min(), KI_max=kv.max(), KI_mean=kv.mean(), diff_KI_minus_WT=kv.mean() - wv.mean(), separation=separation(wv, kv),
                             WT_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(w['animal_id'], wv)), KI_values=';'.join('%s:%.6g' % (i, v) for i, v in zip(ki['animal_id'], kv)))
                    pr = exact_perm(wv, kv)
                    if pr: r.update(pr)
                    else: r['note'] = 'fewer than 2 mice in an arm: permutation not estimable'
                else: r['note'] = 'an arm is empty'
                test_rows.append(r)
    TT = pd.DataFrame(test_rows)
    pcs_t = TT[TT['role'] == 'PC_check'].copy(); pcs_t['PC'] = pcs_t['measure'].str.extract(r'PC(\d+)')[0].astype(int)
    by_pc = pcs_t.pivot_table(index=['population', 'PC'], columns='age_weeks', values='separation', aggfunc='first').reset_index(); by_pc.columns = [str(c) if not isinstance(c, (int, np.integer)) else 'sep_%dW' % c for c in by_pc.columns]
    # ---- write (exclusive create)
    def wr(df, name): df.to_csv(out + name, sep='\t', index=False, float_format='%.10g', mode='x')
    wr(MT, 'mouse_table.tsv'); wr(TT, 'kivswt_within_age.tsv'); wr(by_pc, 'pc_check_by_pc.tsv'); wr(pd.DataFrame(CHECKS, columns=['check_passed', 'items']), 'checks.tsv')
    prim = TT[TT['family'] == 'primary_24']
    nsep = pcs_t.groupby('population')['separation'].apply(lambda s: int(s.isin(['KI_higher', 'KI_lower']).sum())).to_dict()
    res = dict(mode='real' if a.real else 'fixture', populations=POPS, floor_nuclei=FLOOR, mice_included={p: int(MT[(MT.population == p) & MT.included].shape[0]) for p in POPS}, mice_below_floor=MT[~MT['included']][['population', 'animal_id']].values.tolist(),
               checks=len(CHECKS), tests=len(TT), pc_tests_per_population=int(len(pcs_t) // len(POPS)), pc_tests_with_complete_separation=nsep, primary_tests=int(len(prim)), primary_tests_with_complete_separation=int(prim['separation'].isin(['KI_higher', 'KI_lower']).sum()),
               note='3 v 3 complete separation has probability 2/20 = 0.10 per test under the null, so about 10% of the PC tests separate by chance; descriptive, not counted (S11 = C2)',
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, script_sha256=sha256(os.path.abspath(__file__)), seconds=round(time.time() - t0, 2))
    json.dump(res, open(out + 'RESULT.json', 'x'), indent=1, sort_keys=True)
    print(out)


if __name__ == '__main__': main()
#%% debug
