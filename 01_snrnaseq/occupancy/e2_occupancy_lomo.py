"""E2 (pre-specified): Purkinje lane occupancy without double dipping..

Base: stage_a of far_set.py (254eb234), the same lane PCs, lane direction and WT-percentile rule.
  gate  the full-data direction, WT p90 thresholds and per-mouse counts beyond p90 must reproduce the accepted
        far-set run (real_v02 RECORD.json 7ed5c0c1, OCCUPANCY.tsv 2823a594) before anything else is computed.
  (a)   leave-one-mouse-out, every age and mouse: direction from the other 5 mice, WT p90 from the remaining WT
        mice's pooled nuclei, the left-out mouse scored.
  (b)   the six 12-week mice (out-of-sample projection of the frozen discovery PCA, phase4_12w_full) scored on the
        full-data direction and WT p90 of 5 w and of 18 w; basis checked against the stored per-mouse means first.
  (c)   each KI mouse at 18/24/30 w: one Gaussian vs a two-Gaussian mixture with one component fixed to the WT
        mean and SD (ddof 0, pooled WT nuclei) of that age; primary on the full-data scores, LOMO as sensitivity.
Read-only; writes only a NEW real_v01/ next to this script, exclusive create.

  python -I -B 01_snrnaseq/occupancy/e2_occupancy_lomo.py
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
import hashlib
import io
import json
import math

import numpy as np
import pandas as pd

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/e2_occupancy_lomo'))
NOTE = ('Exploratory; E2 occupancy on the Purkinje genotype lanes. Lane PCs were chosen on all 24 discovery mice '
        '(kept fixed, as instructed); the discovery PCA was fitted on all 24 mice.')
AGES = (5, 18, 24, 30)
LANE_CAND = (5, 6, 11, 30)
EXPECTED_LANES = {5: (5, 6, 11), 18: (5, 6, 11, 30), 24: (5, 6, 11, 30), 30: (5, 6, 30)}
P90 = 0.90
W12_MICE = ('2061', '2062', '22018', '2063', '2065', '22019')
SD_FLOOR = 0.05            # mixture: free component SD never below 0.05 x the WT SD
EM_TOL, EM_MAX = 1e-10, 20000
FAR = _release_paths.expand('${REPO_ROOT}/results/far_set_ratios/real_v02/')
W12 = _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/')
PINS = {
    'nucleus_pcs': (_release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'), 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b'),
    'nucleus_meta': (_release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'), '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721'),
    'unit_meta': (_release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'), 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657'),
    'pc_check': (_release_paths.expand('${REPO_ROOT}/results/s12_lanes/pc_check_by_pc.tsv'), '27701bd432255403c50f399ad9a9f98a374b8cdccdfb2e8c22e87849b3bc4b31'),
    'far_record': (FAR + 'RECORD.json', '7ed5c0c184ca977afe71476e6e92afbcd2a75b7e3e1f94536687b64d13df0ea4'),
    'far_occupancy': (FAR + 'OCCUPANCY.tsv', '2823a594148d1b625a4938d3cd8ecb12112a86478fd0b0535b1cf2bcaebf9b1e'),
    'w12_pcs': (W12 + 'Purkinje_cell_PCs_12w.npy', '3c2eb711b5f3c6777c2f7deb43789337633132194bb09d168544d312ddcfb310'),
    'w12_rows': (W12 + 'Purkinje_cell_projection_12w.tsv.gz', '26161c549f37faa58567eca24064603aedec955c5d875996754506e1359d15e4'),
    'w12_means': (W12 + 'Purkinje_cell_PC_means_by_mouse_12w.tsv', '4613c11497fe055297a540595a7b855849ba9595e0ca78deabcb15a4c7ea9158'),
    'ref24_means': (W12 + 'Purkinje_cell_PC_means_by_mouse_24_reference.tsv', '265efe0b5d881c8bddd583d58d64bc7af89e5fb5e9614b603b9fac520ef4e23c'),
    'animals': (_release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'), '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'),
}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def wtsv(df, name):
    buf = io.StringIO()
    buf.write('# ' + NOTE + '\n')
    df.to_csv(buf, sep='\t', index=False, float_format='%.6g')
    with open(OUT / name, 'x') as f:
        f.write(buf.getvalue())


def unit(dvec):
    return dvec / np.linalg.norm(dvec)


def direction(means, ki, wt, lanes):  # unit vector of mean(KI per-mouse means) - mean(WT per-mouse means) on the lanes
    d = np.array([np.mean([means[m][k] for m in ki]) - np.mean([means[m][k] for m in wt]) for k in lanes])
    return unit(d)


def norm_logpdf(x, mu, sd):
    return -0.5 * np.log(2 * np.pi * sd * sd) - (x - mu) ** 2 / (2 * sd * sd)


def one_gaussian(x):
    mu, sd = float(x.mean()), float(x.std(ddof=0))
    ll = float(norm_logpdf(x, mu, sd).sum())
    return dict(mu1=mu, sd1=sd, ll1=ll, bic1=2 * math.log(x.size) - 2 * ll)


def mixture(x, mu0, sd0):  # w * N(mu0, sd0) fixed + (1 - w) * N(mu, sd) free; EM, deterministic starts, best LL kept
    l0 = norm_logpdf(x, mu0, sd0)
    floor = SD_FLOOR * sd0
    best = None
    for w_init in (0.1, 0.5, 0.9):
        for mu_init in (float(np.median(x)), float(np.quantile(x, 0.75))):
            w, mu, sd = w_init, mu_init, max(float(x.std(ddof=0)), floor)
            prev, it, converged = -np.inf, 0, False
            while it < EM_MAX:
                it += 1
                a = np.log(w) + l0
                b = np.log1p(-w) + norm_logpdf(x, mu, sd)
                top = np.maximum(a, b)
                lse = top + np.log(np.exp(a - top) + np.exp(b - top))
                ll = float(lse.sum())
                r0 = np.exp(a - lse)
                r1 = 1.0 - r0
                w = float(np.clip(r0.mean(), 1e-12, 1 - 1e-12))
                if r1.sum() > 0:
                    mu = float((r1 * x).sum() / r1.sum())
                    sd = max(float(math.sqrt((r1 * (x - mu) ** 2).sum() / r1.sum())), floor)
                if abs(ll - prev) < EM_TOL:
                    converged = True
                    break
                prev = ll
            a = np.log(w) + l0
            b = np.log1p(-w) + norm_logpdf(x, mu, sd)
            top = np.maximum(a, b)
            ll = float((top + np.log(np.exp(a - top) + np.exp(b - top))).sum())
            fit = dict(w_WT=w, mu2=mu, sd2=sd, ll2=ll, bic2=3 * math.log(x.size) - 2 * ll, em_iterations=it,
                       converged=converged, start=f'w{w_init}_mu{mu_init:.6g}')
            if best is None or ll > best['ll2'] + 1e-9:
                best = fit
    best['free_on_WT'] = bool(abs(best['mu2'] - mu0) < 0.1 * sd0 and abs(best['sd2'] - sd0) < 0.1 * sd0)
    best['w_at_bound'] = bool(best['w_WT'] < 1e-3 or best['w_WT'] > 1 - 1e-3)
    best['sd2_at_floor'] = bool(best['sd2'] <= floor * (1 + 1e-9))
    return best


def main():
    P = {}
    for key, (rel, want) in PINS.items():
        if sha(ROOT / rel) != want:
            raise ValueError(f'input hash differs: {key} {rel}')
        P[key] = ROOT / rel
    rec = {'note': NOTE, 'script_sha256': sha(Path(__file__)), 'inputs': {rel: h for rel, h in PINS.values()},
           'numpy': np.__version__, 'pandas': pd.__version__}

    # discovery nuclei, as in far_set stage_a
    X = np.load(P['nucleus_pcs'])
    nm = pd.read_csv(P['nucleus_meta'], sep='\t')
    um = pd.read_csv(P['unit_meta'], sep='\t', dtype={'animal_id': str})
    assert X.ndim == 2 and X.shape[0] == len(nm) and X.shape[1] >= 30 and np.isfinite(X).all()
    assert nm.obs_name.is_unique and (um.index.to_numpy() == um.unit_index.to_numpy()).all()
    pk = um[um.cell_type_raw == 'Purkinje cell']
    assert len(pk) == 24 and set(pk.age_weeks) == set(AGES) and (pk.groupby(['age_weeks', 'genotype']).size() == 3).all()
    cell = nm.merge(pk[['unit_index', 'animal_id', 'age_weeks', 'genotype', 'sex']], on='unit_index', how='left', validate='many_to_one')
    assert cell.animal_id.notna().all() and (cell.obs_name.to_numpy() == nm.obs_name.to_numpy()).all()
    mice = pk[['animal_id', 'age_weeks', 'genotype', 'sex']].reset_index(drop=True)
    rows_of = {m: (cell.animal_id == m).to_numpy() for m in mice.animal_id}
    means = {m: {k: float(X[rows_of[m], k - 1].mean()) for k in range(1, 31)} for m in mice.animal_id}

    pcc = pd.read_csv(P['pc_check'], sep='\t')
    pcc = pcc[pcc.population == 'Purkinje_cell'].set_index('PC')
    for a in AGES:
        lanes = tuple(k for k in LANE_CAND if pcc.loc[k, f'sep_{a}W'] != 'none')
        assert lanes == EXPECTED_LANES[a], ('lane PCs differ from the pre-registered set', a, lanes)

    def scores(rows, a, u):
        return X[rows][:, [k - 1 for k in EXPECTED_LANES[a]]] @ u

    # full-data direction and thresholds; gate against the accepted far-set run
    far = json.loads(P['far_record'].read_text())
    occ_old = pd.read_csv(P['far_occupancy'], sep='\t', comment='#', dtype={'mouse': str})
    full_u, full_thr, gate_dev = {}, {}, 0.0
    for a in AGES:
        am = mice[mice.age_weeks == a]
        ki, wt = am[am.genotype == 'KI'].animal_id.tolist(), am[am.genotype == 'WT'].animal_id.tolist()
        u = direction(means, ki, wt, EXPECTED_LANES[a])
        full_u[a] = u
        wsc = np.concatenate([scores(rows_of[m], a, u) for m in wt])
        full_thr[a] = float(np.quantile(wsc, P90))
        old_u = np.array([far['lane_unit_vector'][str(a)][f'PC{k}'] for k in EXPECTED_LANES[a]])
        gate_dev = max(gate_dev, float(np.abs(old_u - u).max()), abs(far['wt_quantile_thresholds'][str(a)]['0.9'] - full_thr[a]))
        for m in ki + wt:
            n_new = int((scores(rows_of[m], a, u) > full_thr[a]).sum())
            n_old = int(occ_old.loc[(occ_old.mouse == m) & (occ_old.age == a), 'n_beyond_WT_p90'].iloc[0])
            if n_new != n_old:
                raise ValueError(f'gate: n beyond p90 differs from the accepted far-set run: {a} w {m} {n_new} vs {n_old}')
    if gate_dev > 1e-12:
        raise ValueError(f'gate: direction or p90 differs from the accepted far-set run by {gate_dev}')
    rec['gate_max_abs_deviation_vs_far_set_real_v02'] = gate_dev

    # (a) leave one mouse out
    lomo = []
    for a in AGES:
        am = mice[mice.age_weeks == a]
        for _, m in am.iterrows():
            rest = am[am.animal_id != m.animal_id]
            ki, wt = rest[rest.genotype == 'KI'].animal_id.tolist(), rest[rest.genotype == 'WT'].animal_id.tolist()
            u = direction(means, ki, wt, EXPECTED_LANES[a])
            wsc = np.concatenate([scores(rows_of[w], a, u) for w in wt])
            thr = float(np.quantile(wsc, P90))
            s = scores(rows_of[m.animal_id], a, u)
            s_full = scores(rows_of[m.animal_id], a, full_u[a])
            lomo.append(dict(age=a, mouse=m.animal_id, genotype=m.genotype, sex=m.sex, n_nuclei=s.size,
                             n_KI_mice_in_direction=len(ki), n_WT_mice_in_direction=len(wt), n_WT_nuclei_in_threshold=wsc.size,
                             lomo_p90=thr, lomo_n_beyond=int((s > thr).sum()), lomo_share=float((s > thr).mean()),
                             full_p90=full_thr[a], full_n_beyond=int((s_full > full_thr[a]).sum()),
                             full_share=float((s_full > full_thr[a]).mean()),
                             cos_lomo_vs_full=float(u @ full_u[a]),
                             **{f'lomo_u_PC{k}': float(v) for k, v in zip(EXPECTED_LANES[a], u)}))
    lomo = pd.DataFrame(lomo)

    # (b) 12-week mice; basis checks first
    ref = pd.read_csv(P['ref24_means'], sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    basis_dev = max(abs(means[m][k] - float(ref.loc[m, f'PC{k}'])) for m in mice.animal_id for k in range(1, 31))
    if basis_dev > 1e-9:
        raise ValueError(f'basis check: discovery per-mouse means differ from the 24-mouse reference by {basis_dev}')
    X12 = np.load(P['w12_pcs'])
    r12 = pd.read_csv(P['w12_rows'], sep='\t', dtype={'animal_id': str})
    assert X12.shape == (len(r12), 30) and np.isfinite(X12).all() and set(r12.animal_id) == set(W12_MICE)
    m12 = pd.read_csv(P['w12_means'], sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    w12_dev = max(abs(float(X12[(r12.animal_id == m).to_numpy(), k - 1].mean()) - float(m12.loc[m, f'PC{k}']))
                  for m in W12_MICE for k in range(1, 31))
    if w12_dev > 1e-9:
        raise ValueError(f'basis check: 12-week per-mouse means differ from the stored table by {w12_dev}')
    animals = pd.read_csv(P['animals'], sep='\t', dtype={'animal_id': str}).set_index('animal_id')
    for m in W12_MICE:
        g = animals.loc[m, 'genotype']
        assert int(animals.loc[m, 'age_weeks']) == 12 and set(r12.loc[r12.animal_id == m, 'genotype']) == {g}
        assert int((r12.animal_id == m).sum()) == int(m12.loc[m, 'n_nuclei'])
    rec['basis_check_max_abs_deviation'] = {'discovery_24_reference': basis_dev, 'w12_stored_means': w12_dev}
    w12 = []
    for ref_age in (5, 18):
        u, thr, lanes = full_u[ref_age], full_thr[ref_age], EXPECTED_LANES[ref_age]
        for m in W12_MICE:
            s = X12[(r12.animal_id == m).to_numpy()][:, [k - 1 for k in lanes]] @ u
            w12.append(dict(reference_age=ref_age, mouse=m, genotype=animals.loc[m, 'genotype'], sex=animals.loc[m, 'sex'],
                            n_nuclei=s.size, ref_p90=thr, n_beyond=int((s > thr).sum()), share=float((s > thr).mean()),
                            lane_score_median=float(np.median(s))))
    w12 = pd.DataFrame(w12)

    # (c) one Gaussian vs WT-anchored mixture, KI mice at 18/24/30 w
    mix = []
    for a in (18, 24, 30):
        am = mice[mice.age_weeks == a]
        for _, m in am[am.genotype == 'KI'].iterrows():
            for version in ('full', 'lomo'):
                if version == 'full':
                    u, wt = full_u[a], am[am.genotype == 'WT'].animal_id.tolist()
                else:
                    rest = am[am.animal_id != m.animal_id]
                    wt = rest[rest.genotype == 'WT'].animal_id.tolist()
                    u = direction(means, rest[rest.genotype == 'KI'].animal_id.tolist(), wt, EXPECTED_LANES[a])
                wsc = np.concatenate([scores(rows_of[w], a, u) for w in wt])
                mu0, sd0 = float(wsc.mean()), float(wsc.std(ddof=0))
                x = scores(rows_of[m.animal_id], a, u)
                g1, g2 = one_gaussian(x), mixture(x, mu0, sd0)
                mix.append(dict(version=version, age=a, mouse=m.animal_id, sex=m.sex, n_nuclei=x.size,
                                WT_mean=mu0, WT_sd=sd0, n_WT_nuclei=wsc.size, **g1, **g2,
                                delta_bic_2_minus_1=g2['bic2'] - g1['bic1']))
    mix = pd.DataFrame(mix)

    OUT.mkdir(parents=False, exist_ok=False)
    wtsv(lomo, 'E2a_LOMO.tsv')
    wtsv(w12, 'E2b_W12_SHARES.tsv')
    wtsv(mix, 'E2c_MIXTURE.tsv')
    rec.update(lane_pcs={str(a): list(EXPECTED_LANES[a]) for a in AGES},
               full_unit_vector={str(a): {f'PC{k}': float(v) for k, v in zip(EXPECTED_LANES[a], full_u[a])} for a in AGES},
               full_wt_p90={str(a): full_thr[a] for a in AGES}, quantile='numpy default (linear), pooled WT nuclei; beyond = strictly greater',
               mixture={'free_sd_floor_x_WT_sd': SD_FLOOR, 'em_tol': EM_TOL, 'em_max': EM_MAX,
                        'starts': 'w in {0.1, 0.5, 0.9} x mu in {median, q75}; best log-likelihood kept',
                        'bic': 'k ln n - 2 LL; k = 2 (one Gaussian), 3 (mixture)', 'WT_sd': 'ddof 0'})
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print(json.dumps({'gate_dev': gate_dev, 'basis_dev': basis_dev, 'w12_dev': w12_dev,
                      **{p.name: sha(p) for p in files}}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
