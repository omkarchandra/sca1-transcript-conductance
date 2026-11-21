
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Phase 4 A3: bimodality check of the KI spread (pre-specified). CODE ONLY. NOT APPROVED TO RUN: it answers an open decision of <AUTHOR>'s. DESCRIPTIVE counts only; no p-values; no biological reading.
# HARD GATE: every mode (--check, --full, --replay) refuses to open any real file unless --go-file PATH names a JSON that (a) holds <AUTHOR>'s written yes (non-empty 'author_written_yes'), (b) is bound to THESE bytes ('script_sha256' == sha256 of this file) and
# (c) says who relayed it ('relayed_by' and 'relayed_at_cdt'). Without it the script stops before reading anything.
# Data (once released): per-nucleus lane scores of the frozen 24-mouse S07 projection (nucleus_PCs.npy via nucleus_metadata.unit_index and the stage A unit table), at 24 W and 30 W only, on the 12 pre-declared lanes (Purkinje 5, 6, 11, 30; Bergmann 2, 4, 5, 6, 8, 9, 20, 29;
# lane k = column k-1), analysed PER MOUSE (3 KI + 3 WT per age). A mouse with fewer than 20 nuclei is flagged and EXCLUDED (no statistics). Labels and ages come from the pinned animal manifest; the stage A unit table must agree.
# Per mouse x lane: (a) 1- vs 2-component Gaussian mixture, delta_BIC = BIC(1) - BIC(2), sklearn.mixture.GaussianMixture(covariance_type='full', n_init=10, random_state=SEED, max_iter=500; tol and reg_covar at their defaults), the same fixed seed for every fit;
#                   (b) Sarle's bimodality coefficient BC = (G1^2 + 1) / (G2 + 3(n-1)^2/((n-2)(n-3))) with the unbiased skewness G1 and unbiased excess kurtosis G2 (descriptive; 5/9 is the usual reference value, shown as a plain flag).
# 'TWO-GROUP FAVOURED' call: delta_BIC > 10 (strict) AND smaller mixture weight >= 0.10 AND the two component means >= 2 pooled within-component SDs apart; otherwise 'one group'.
#   Pooled within-component SD = sqrt(w1*var1 + w2*var2) (weighted by the mixture weights). This definition is MINE where the rule is silent; the other reading, sqrt((var1 + var2)/2), is a question for review.
# Lane-age call 'yes' (= KI-specific two-group structure) only if >= 2 of 3 KI mice are flagged AND 0 of 3 WT mice are flagged at that lane and age. If an excluded mouse leaves fewer than 3 KI or 3 WT mice analysed, the lane-age is 'NOT_EVALUABLE' (counted separately, never as yes).
# Report: the number of lane-ages meeting the rule out of 24 (12 lanes x 2 ages). Falsifier: none meets it.
# Run (only after the go-file exists): python -I -B -I ${REPO_ROOT}/01_snrnaseq/spread/phase4_a3_bimodality.py --go-file PATH --check | --full | --replay DIR      (python -I -B; exclusive-create outputs; no deletion)
import argparse, csv, hashlib, json, math, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd
import scipy, sklearn
from scipy.stats import kurtosis, skew
from sklearn.mixture import GaussianMixture

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
TARGETS = {'Purkinje_cell': (5, 6, 11, 30), 'Bergmann_glia': (2, 4, 5, 6, 8, 9, 20, 29)}
CELLTYPE = {'Purkinje_cell': 'Purkinje cell', 'Bergmann_glia': 'Bergmann glia'}
AGES = (24, 30); FLOOR = 20
DELTA_BIC = 10.0; MIN_WEIGHT = 0.10; MIN_SEP = 2.0; BC_REFERENCE = 5 / 9
SEED = 20261007; N_INIT = 10; MAX_ITER = 500
CFG = dict(units=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'), core=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01'), animal_manifest=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'))
PINS = {  # relative path -> full sha256 (the same pins the accepted A1 script used for the same files)
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'): 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'): 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_PCs.npy'): 'ed535ac9cdb36c51498155641172de92d9b6daa65e611ffa81a705b1a573ea13',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): '23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540',
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'}
yn = lambda b: 'yes' if b else 'no'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def check_go(go_file, script_path=None):
    """The gate: <AUTHOR>'s written yes, bound to these bytes, with who relayed it and when. Raises before any real file is opened."""
    g = json.loads(Path(go_file).read_text()); me = sha(script_path or Path(__file__))
    assert isinstance(g.get('author_written_yes'), str) and g['author_written_yes'].strip(), 'no written yes from <AUTHOR> in the go-file'
    assert g.get('script_sha256') == me, ('the go-file is bound to other bytes than this script', g.get('script_sha256'), me)
    assert str(g.get('relayed_by', '')).strip() and str(g.get('relayed_at_cdt', '')).strip(), 'the go-file must say who relayed the yes and when'
    return g


def check_pins(pins, root):
    got = {}
    for rel, h in pins.items():
        got[rel] = sha(Path(root) / rel); assert got[rel] == h, ('pinned file changed', rel)
    return got


def sarle_bc(x):
    """Sarle's bimodality coefficient with the unbiased skewness and unbiased excess kurtosis (descriptive)."""
    n = len(x); g1 = float(skew(x, bias=False)); k = float(kurtosis(x, fisher=True, bias=False))
    return (g1 ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def mixture_stats(x):
    """1- vs 2-component Gaussian mixture on one mouse x lane sample. Component 'small' = the smaller mixture weight."""
    X = np.asarray(x, float).reshape(-1, 1)
    g1 = GaussianMixture(n_components=1, covariance_type='full', n_init=N_INIT, random_state=SEED, max_iter=MAX_ITER).fit(X)
    g2 = GaussianMixture(n_components=2, covariance_type='full', n_init=N_INIT, random_state=SEED, max_iter=MAX_ITER).fit(X)
    w = g2.weights_.ravel(); mu = g2.means_.ravel(); var = g2.covariances_.ravel(); i = np.argsort(w, kind='stable'); a, b = int(i[0]), int(i[1])
    pooled = float(np.sqrt((w * var).sum()))
    return dict(bic_1=float(g1.bic(X)), bic_2=float(g2.bic(X)), delta_bic=float(g1.bic(X) - g2.bic(X)), w_small=float(w[a]), mu_small=float(mu[a]), mu_large=float(mu[b]), sd_small=float(np.sqrt(var[a])), sd_large=float(np.sqrt(var[b])),
                pooled_sd=pooled, separation_sd=float(abs(mu[a] - mu[b]) / pooled), converged_1=bool(g1.converged_), converged_2=bool(g2.converged_))


def two_group_call(delta_bic, w_small, sep):
    return bool(delta_bic > DELTA_BIC and w_small >= MIN_WEIGHT and sep >= MIN_SEP)


def lane_age_call(ki_flags, wt_flags):
    if len(ki_flags) != 3 or len(wt_flags) != 3: return 'NOT_EVALUABLE'
    return 'yes' if sum(ki_flags) >= 2 and sum(wt_flags) == 0 else 'no'


def load_discovery(cfg, pop, manifest):
    """Nucleus lane scores of one population with the mouse, genotype and age of every row (found by unit_index; labels and ages from the manifest, cross-checked against the stage A unit table)."""
    name = f'discovery__{pop}'; nm = pd.read_csv(Path(cfg['core']) / name / 'nucleus_metadata.tsv.gz', sep='\t'); pcs = np.load(Path(cfg['core']) / name / 'nucleus_PCs.npy')
    units = pd.read_csv(cfg['units'], sep='\t', dtype={'animal_id': str}); assert (units.unit_index.to_numpy() == np.arange(len(units))).all()
    u = units.iloc[nm.unit_index.to_numpy()]; g = manifest.set_index('animal_id'); assert len(nm) == len(pcs) and (u.cell_type_raw == CELLTYPE[pop]).all(), pop
    assert (u.genotype.to_numpy() == g.genotype.loc[u.animal_id].to_numpy()).all() and (u.age_weeks.to_numpy() == g.age_weeks.loc[u.animal_id].to_numpy()).all(), ('stage A unit table disagrees with the animal manifest', pop)
    return pcs, u.animal_id.to_numpy().astype(object), g.genotype.loc[u.animal_id].to_numpy().astype(object), g.age_weeks.loc[u.animal_id].to_numpy().astype(int)


def build(cfg):
    manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}); rows = []
    for pop, ks in TARGETS.items():
        pcs, ids, geno, age = load_discovery(cfg, pop, manifest)
        for k in ks:
            for a in AGES:
                mice = sorted({(m, g_) for m, g_, ag in zip(ids, geno, age) if ag == a}); assert sorted(g_ for _m, g_ in mice) == ['KI'] * 3 + ['WT'] * 3, ('need 3 KI and 3 WT mice at this age', pop, a)
                for m, g_ in mice:
                    x = pcs[(ids == m) & (age == a), k - 1]; n = int(len(x)); base = dict(population=pop, PC=k, age_weeks=a, animal_id=m, genotype=g_, n_nuclei=n)
                    if n < FLOOR:
                        rows.append(dict(base, analysed='no', two_group_call='EXCLUDED_LT20')); continue
                    s = mixture_stats(x); bc = sarle_bc(x); c1, c2, c3 = s['delta_bic'] > DELTA_BIC, s['w_small'] >= MIN_WEIGHT, s['separation_sd'] >= MIN_SEP
                    rows.append(dict(base, analysed='yes', **s, sarle_bc=bc, bc_gt_5_9_descriptive=yn(bc > BC_REFERENCE), crit_delta_bic_gt_10=yn(c1), crit_weight_ge_0_10=yn(c2), crit_separation_ge_2_sd=yn(c3), two_group_call='two_group' if two_group_call(s['delta_bic'], s['w_small'], s['separation_sd']) else 'one_group'))
    mouse = pd.DataFrame(rows); la = []
    for (pop, k, a), grp in mouse.groupby(['population', 'PC', 'age_weeks'], sort=False):
        ki = grp[(grp.genotype == 'KI') & (grp.analysed == 'yes')]; wt = grp[(grp.genotype == 'WT') & (grp.analysed == 'yes')]
        call = lane_age_call([c == 'two_group' for c in ki.two_group_call], [c == 'two_group' for c in wt.two_group_call])
        la.append(dict(population=pop, PC=int(k), age_weeks=int(a), n_KI_analysed=len(ki), n_WT_analysed=len(wt), KI_flagged=int((ki.two_group_call == 'two_group').sum()), WT_flagged=int((wt.two_group_call == 'two_group').sum()), KI_specific_two_group=call))
    lane_age = pd.DataFrame(la)
    return mouse, lane_age, summarize(mouse, lane_age)


def summarize(mouse, lane_age):
    out = {}
    for name, m, l in [(p, mouse[mouse.population == p], lane_age[lane_age.population == p]) for p in TARGETS] + [('both', mouse, lane_age)]:
        met = int((l.KI_specific_two_group == 'yes').sum()); ne = int((l.KI_specific_two_group == 'NOT_EVALUABLE').sum())
        calls = {f'{g}_{a}W': dict(analysed=int(((m.genotype == g) & (m.age_weeks == a) & (m.analysed == 'yes')).sum()), two_group=int(((m.genotype == g) & (m.age_weeks == a) & (m.two_group_call == 'two_group')).sum())) for g in ('KI', 'WT') for a in AGES}
        out[name] = dict(n_lane_ages=int(len(l)), n_meeting_rule=met, n_not_evaluable=ne, n_excluded_mouse_lane_rows=int((m.two_group_call == 'EXCLUDED_LT20').sum()), none_meet_rule=bool(met == 0), two_group_calls_by_genotype_and_age=calls,
                         meeting_by_age={f'{a}W': int(((l.age_weeks == a) & (l.KI_specific_two_group == 'yes')).sum()) for a in AGES})
    return out


def write_outputs(out_data, out_docs, mouse, lane_age, summary, pins, seconds):
    out_data.mkdir(parents=True); out_docs.mkdir(parents=True)
    mouse.to_csv(out_data / 'PER_MOUSE_LANE.tsv', sep='\t', index=False); lane_age.to_csv(out_data / 'LANE_AGE_CALLS.tsv', sep='\t', index=False); (out_data / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out_data / 'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(out_data.iterdir()) if p.is_file()))
    (out_docs / 'RECORD.json').write_text(json.dumps(dict(seconds=seconds, summary=summary, pins=pins, script_sha256=sha(Path(__file__)), versions=dict(sklearn=sklearn.__version__, scipy=scipy.__version__, numpy=np.__version__, pandas=pd.__version__),
                                                         settings=dict(seed=SEED, n_init=N_INIT, max_iter=MAX_ITER, covariance_type='full', delta_bic_strict_greater_than=DELTA_BIC, min_weight=MIN_WEIGHT, min_separation_sd=MIN_SEP, floor_nuclei=FLOOR, ages=list(AGES)),
                                                         note='descriptive counts only; rule = pre-registered planning note section 18 (proposed); pooled within-component SD = sqrt(w1*var1 + w2*var2)'), indent=2) + '\n')


def replay(out_data, cfg=CFG):
    """Fresh-process replay. (1) Every saved statistic is re-derived in PLAIN PYTHON (math only) from the raw nucleus scores and the SAVED mixture parameters: n, BIC(1) in closed form, BIC(2) from the log-likelihood at the saved parameters, delta_BIC, pooled SD,
    separation, Sarle's BC with hand-written moments, every criterion and call, and the lane-age tallies and summary counts from the saved table. (2) The mixture fits are repeated with the same seed and must reproduce the saved parameters to 1e-8 (determinism)."""
    out_data = Path(out_data); mouse = pd.read_csv(out_data / 'PER_MOUSE_LANE.tsv', sep='\t', dtype={'animal_id': str}, keep_default_na=False); bad = []; manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}); cache = {}
    assert len(mouse) == sum(len(v) for v in TARGETS.values()) * 2 * 6
    def logn(x, mu, var): return -0.5 * (math.log(2 * math.pi * var) + (x - mu) ** 2 / var)
    for pop in TARGETS: cache[pop] = load_discovery(cfg, pop, manifest)
    for r in mouse.itertuples():
        pcs, ids, geno, age = cache[r.population]; x = [float(v) for v in pcs[(ids == r.animal_id) & (age == r.age_weeks), r.PC - 1]]; n = len(x); ok = n == r.n_nuclei and r.genotype == manifest.set_index('animal_id').genotype[r.animal_id]
        if n < FLOOR: ok = ok and r.analysed == 'no' and r.two_group_call == 'EXCLUDED_LT20'
        else:
            m = sum(x) / n; v1 = sum((t - m) ** 2 for t in x) / n + 1e-6; bic1 = -2 * sum(logn(t, m, v1) for t in x) + 2 * math.log(n)
            ws, wl = float(r.w_small), 1 - float(r.w_small); ll = 0.0
            for t in x: ll += math.log(ws * math.exp(logn(t, float(r.mu_small), float(r.sd_small) ** 2)) + wl * math.exp(logn(t, float(r.mu_large), float(r.sd_large) ** 2)))
            bic2 = -2 * ll + 5 * math.log(n); pooled = math.sqrt(ws * float(r.sd_small) ** 2 + wl * float(r.sd_large) ** 2); sep = abs(float(r.mu_small) - float(r.mu_large)) / pooled; s = math.sqrt(sum((t - m) ** 2 for t in x) / (n - 1))
            g1 = n / ((n - 1) * (n - 2)) * sum(((t - m) / s) ** 3 for t in x); g2 = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * sum(((t - m) / s) ** 4 for t in x) - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)); bc = (g1 ** 2 + 1) / (g2 + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))
            rel = lambda a, b: abs(a - b) <= 1e-6 * max(1.0, abs(b)); call = 'two_group' if (bic1 - bic2 > DELTA_BIC and ws >= MIN_WEIGHT and sep >= MIN_SEP) else 'one_group'
            ok = ok and r.analysed == 'yes' and rel(bic1, float(r.bic_1)) and rel(bic2, float(r.bic_2)) and rel(bic1 - bic2, float(r.delta_bic)) and rel(pooled, float(r.pooled_sd)) and rel(sep, float(r.separation_sd)) and rel(bc, float(r.sarle_bc)) and r.two_group_call == call \
                and r.bc_gt_5_9_descriptive == yn(bc > BC_REFERENCE) and r.crit_delta_bic_gt_10 == yn(bic1 - bic2 > DELTA_BIC) and r.crit_weight_ge_0_10 == yn(ws >= MIN_WEIGHT) and r.crit_separation_ge_2_sd == yn(sep >= MIN_SEP)
            fit = mixture_stats(x); ok = ok and all(abs(fit[k_] - float(getattr(r, k_))) <= 1e-8 * max(1.0, abs(fit[k_])) for k_ in ('bic_1', 'bic_2', 'w_small', 'mu_small', 'mu_large', 'sd_small', 'sd_large'))
        if not ok: bad.append((r.population, int(r.PC), int(r.age_weeks), r.animal_id))
    la_saved = pd.read_csv(out_data / 'LANE_AGE_CALLS.tsv', sep='\t', keep_default_na=False); bad_la = []
    for r in la_saved.itertuples():
        sub = mouse[(mouse.population == r.population) & (mouse.PC == r.PC) & (mouse.age_weeks == r.age_weeks) & (mouse.analysed == 'yes')]; ki = [c == 'two_group' for c in sub[sub.genotype == 'KI'].two_group_call]; wt = [c == 'two_group' for c in sub[sub.genotype == 'WT'].two_group_call]
        want = 'NOT_EVALUABLE' if len(ki) != 3 or len(wt) != 3 else ('yes' if sum(ki) >= 2 and sum(wt) == 0 else 'no')
        if not (r.KI_specific_two_group == want and r.n_KI_analysed == len(ki) and r.n_WT_analysed == len(wt) and r.KI_flagged == sum(ki) and r.WT_flagged == sum(wt)): bad_la.append((r.population, int(r.PC), int(r.age_weeks)))
    summ = json.loads((out_data / 'SUMMARY.json').read_text()); sums_ok = True
    for scope, pops in [(p, [p]) for p in TARGETS] + [('both', list(TARGETS))]:
        l = la_saved[la_saved.population.isin(pops)]; e = summ[scope]
        sums_ok = sums_ok and e['n_lane_ages'] == len(l) and e['n_meeting_rule'] == int((l.KI_specific_two_group == 'yes').sum()) and e['n_not_evaluable'] == int((l.KI_specific_two_group == 'NOT_EVALUABLE').sum()) and e['none_meet_rule'] == (int((l.KI_specific_two_group == 'yes').sum()) == 0)
    res = dict(n_mouse_lane_rows_checked=len(mouse), mouse_lane_rows_with_differences=[list(b) for b in bad], lane_age_rows_with_differences=[list(b) for b in bad_la], summary_counts_equal=bool(sums_ok)); res['all_pass'] = bool(not bad and not bad_la and sums_ok)
    return res


def guarded(out_docs, fn):
    """Run fn; on ANY failure after the output dirs exist, exclusive-create FAILED.json (traceback) in out_docs, then re-raise (no retry, nothing deleted)."""
    try: return fn()
    except BaseException:
        import traceback
        if Path(out_docs).exists():
            with open(Path(out_docs) / 'FAILED.json', 'x') as f: f.write(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), traceback=traceback.format_exc()), indent=2) + '\n')
        raise


def structure_check(cfg):
    """--check (also gated): every population loads, the row mapping asserts hold, every lane column exists, each analysed age has 3 KI + 3 WT mice. No mixture is fitted."""
    manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str})
    for pop, ks in TARGETS.items():
        pcs, ids, geno, age = load_discovery(cfg, pop, manifest); assert pcs.shape[1] >= max(ks)
        for a in AGES: assert sorted(g for _m, g in {(m, g) for m, g, ag in zip(ids, geno, age) if ag == a}) == ['KI'] * 3 + ['WT'] * 3, (pop, a)
    return True


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--go-file', type=Path, required=True); g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--check', action='store_true'); g.add_argument('--full', action='store_true'); g.add_argument('--replay', type=Path)
    a = ap.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'): sys.exit('refusing: the checks are asserts; run without -O')
    go = check_go(a.go_file)  # BEFORE any real file is opened
    pins = check_pins(PINS, ROOT)
    if a.check: print(json.dumps(dict(go_file_ok=True, pins_ok=True, n_pins=len(pins), structure_ok=structure_check(CFG)))); return
    if a.replay:
        res = replay(a.replay); docs = ROOT / 'docs' / a.replay.name
        with open(docs / 'REPLAY.json', 'x') as f: f.write(json.dumps(res, indent=2) + '\n')
        print(json.dumps(res)); sys.exit(0 if res['all_pass'] else 1)
    out_data = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_a3_bimodality'); out_docs = ROOT / _release_paths.expand('${DATA_ROOT}/phase4_a3_bimodality')
    assert not out_data.exists() and not out_docs.exists(), 'outputs exist; never overwrite'
    t0 = time.time(); mouse, lane_age, summary = build(CFG)
    def body():
        write_outputs(out_data, out_docs, mouse, lane_age, summary, dict(pins, go_file=str(a.go_file), go_file_sha256=sha(a.go_file)), time.time() - t0); print(json.dumps(summary))
    guarded(out_docs, body)


if __name__ == '__main__':
    main()

#%% debug
