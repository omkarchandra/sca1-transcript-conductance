
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Phase 4, 12 W held-out check (pre-specified). CODE ONLY until review clears it by hash.
# For the 12 pre-declared PCs (Purkinje 5, 6, 11, 30; Bergmann 2, 4, 5, 6, 8, 9, 20, 29): predicted sign = majority sign of the discovery KI-minus-WT differences at the 4 discovery ages (T3 table, read-only);
# a 2-2 split = no prediction (excluded, listed). 12 W: KI-minus-WT difference of the 3 v 3 per-mouse PC means (from the accepted full run), sign match, complete separation in the predicted direction (strict).
# Summary: k/n sign matches (chance n/2), m/n complete separations (chance 0.05 each: 3 v 3 has 20 labelings, one of them is the predicted direction). Descriptive; NO verdict, no p-value, no biological reading.
# Labels (WT/KI) come from the pinned animal manifest only; the genotype column of the PC-means tables must agree with it. Run with: python -I -B ${REPO_ROOT}/01_snrnaseq/heldout_12w/phase4_12w_heldout_check.py --check | --full | --replay DIR
import argparse, csv, hashlib, json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
SIX = ('2061', '2062', '22018', '2063', '2065', '22019')            # the six 12 W mice (3 WT, 3 KI; labels from the manifest)
TARGETS = {'Purkinje_cell': (5, 6, 11, 30), 'Bergmann_glia': (2, 4, 5, 6, 8, 9, 20, 29)}
AGES = (5, 18, 24, 30)
CFG = dict(t3=ROOT / _release_paths.expand('${REPO_ROOT}/results/s12_lanes/kivswt_within_age.tsv'), pcmeans=ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full'), animal_manifest=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'))
PINS = {  # relative path -> full sha256
    _release_paths.expand('${REPO_ROOT}/results/s12_lanes/kivswt_within_age.tsv'): '01f29ce2b9c9e648b37e478ed408e9a3d276df5240d776e6996fa1cd24565914',
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Purkinje_cell_PC_means_by_mouse_12w.tsv'): '4613c11497fe055297a540595a7b855849ba9595e0ca78deabcb15a4c7ea9158',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Bergmann_glia_PC_means_by_mouse_12w.tsv'): 'f71c72f1408ae0865c019e63e19b845cc5f8616e0bfbee44efdb77fee0fe7f08',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/SHA256SUMS.txt'): '187e91a7ddaac70c026ee90481cbc02cad6a672be1bed843a01a04c5f9b423fd',
    _release_paths.expand('${DATA_ROOT}/phase4_12w_full/RECORD.json'): '6b0e3c442b1e6ed41ecd2d67ee1d49e59149c22cfcafedadaab7956f0b784123'}
SYM = {1: '+', -1: '-', 0: '0'}; PSYM = {1: '+', -1: '-', 0: 'none'}; OBS = {1: 'KI_higher', -1: 'KI_lower', 0: 'none'}
yn = lambda b: 'yes' if b else 'no'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(16 << 20), b''): h.update(b)
    return h.hexdigest()


def check_pins(pins, root):
    got = {}
    for rel, h in pins.items():
        got[rel] = sha(Path(root) / rel); assert got[rel] == h, ('pinned file changed', rel)
    return got


def predicted_sign(diffs):
    """Majority sign of the discovery KI-minus-WT differences; equal numbers of + and - (a 2-2 split) = 0 = no prediction."""
    pos = sum(1 for d in diffs if d > 0); neg = sum(1 for d in diffs if d < 0)
    return 1 if pos > neg else (-1 if neg > pos else 0)


def discovery_diffs(t3, pop, k):
    """{age: KI-minus-WT} of discovery PC k, found by population/measure/age_weeks (no row order assumed); each of the 4 ages exactly once."""
    r = t3[(t3.population == pop) & (t3.measure == f'discovery_PC{k}_mean')]
    assert sorted(r.age_weeks.tolist()) == list(AGES) and (r.role == 'PC_check').all(), (pop, k)
    return {int(a): float(d) for a, d in zip(r.age_weeks, r.diff_KI_minus_WT)}


def score_pc(pred, wt, ki):
    """wt, ki: the three per-mouse PC means of each genotype. diff = mean(KI) - mean(WT). Complete separation is strict (a tie is not separation)."""
    diff = float(np.mean(ki) - np.mean(wt))
    sign = 1 if diff > 0 else (-1 if diff < 0 else 0)
    observed = 1 if min(ki) > max(wt) else (-1 if max(ki) < min(wt) else 0)
    if pred == 0: return dict(diff=diff, sign=sign, observed=observed, sign_match='NA', complete_separation='NA')
    return dict(diff=diff, sign=sign, observed=observed, sign_match=yn(sign == pred), complete_separation=yn(observed == pred))


def groups(manifest, means_by_pop):
    """WT and KI mice from the pinned animal manifest (NOT from the PC-means tables, whose genotype column must agree)."""
    g = manifest.set_index('animal_id').genotype
    wt = [a for a in SIX if g[a] == 'WT']; ki = [a for a in SIX if g[a] == 'KI']
    assert len(wt) == 3 and len(ki) == 3
    for pop, m in means_by_pop.items():
        assert sorted(m.animal_id.tolist()) == sorted(SIX), pop
        assert list(m.set_index('animal_id').genotype.loc[list(SIX)]) == [g[a] for a in SIX], ('PC-means genotype column disagrees with the animal manifest', pop)
    return wt, ki


def build(cfg):
    t3 = pd.read_csv(cfg['t3'], sep='\t'); manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}); rows = []
    means_by_pop = {pop: pd.read_csv(Path(cfg['pcmeans']) / f'{pop}_PC_means_by_mouse_12w.tsv', sep='\t', dtype={'animal_id': str}) for pop in TARGETS}
    wt, ki = groups(manifest, means_by_pop)
    for pop, pcs in TARGETS.items():
        means = means_by_pop[pop].set_index('animal_id')
        for k in pcs:
            d = discovery_diffs(t3, pop, k); pred = predicted_sign([d[a] for a in AGES])
            col = means[f'PC{k}']; w = col.loc[wt].to_numpy(float); q = col.loc[ki].to_numpy(float); r = score_pc(pred, w, q)
            row = dict(population=pop, PC=k); row.update({f'disc_diff_{a}W': d[a] for a in AGES}); row.update({f'disc_sign_{a}W': SYM[int(np.sign(d[a]))] for a in AGES})
            row.update(n_plus=sum(d[a] > 0 for a in AGES), n_minus=sum(d[a] < 0 for a in AGES), predicted_sign=PSYM[pred], excluded=pred == 0,
                       w12_WT_means=';'.join(f'{a}:{v:.6g}' for a, v in zip(wt, w)), w12_KI_means=';'.join(f'{a}:{v:.6g}' for a, v in zip(ki, q)), w12_diff_KI_minus_WT=r['diff'], w12_sign=SYM[r['sign']],
                       observed_complete_separation=OBS[r['observed']], sign_match=r['sign_match'], complete_separation_in_predicted_direction=r['complete_separation'])
            rows.append(row)
    table = pd.DataFrame(rows)
    return table, summarize(table), {pop: means_by_pop[pop].loc[means_by_pop[pop].below_20_nuclei_flag.astype(str).str.lower().eq('true'), 'animal_id'].tolist() for pop in TARGETS}


def summarize(table):
    out = {}
    for name, sub in [(p, table[table.population == p]) for p in TARGETS] + [('both', table)]:
        pred = sub[~sub.excluded]; n = int(len(pred))
        out[name] = dict(n_listed=int(len(sub)), excluded_2_2_split=[f'{p}:PC{k}' for p, k in zip(sub[sub.excluded].population, sub[sub.excluded].PC)], n_with_prediction=n,
                         sign_matches=int((pred.sign_match == 'yes').sum()), chance_expected_sign_matches=n / 2,
                         complete_separations_in_predicted_direction=int((pred.complete_separation_in_predicted_direction == 'yes').sum()), chance_expected_complete_separations=n * 0.05)
    return out


def write_outputs(out_data, out_docs, table, summary, below20, pins, seconds):
    out_data.mkdir(parents=True); out_docs.mkdir(parents=True)
    table.to_csv(out_data / 'HELDOUT_BY_PC.tsv', sep='\t', index=False); (out_data / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out_data / 'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(out_data.iterdir()) if p.is_file()))
    (out_docs / 'RECORD.json').write_text(json.dumps(dict(seconds=seconds, summary=summary, mice_below_20_nuclei=below20, pins=pins, script_sha256=sha(Path(__file__)),
                                                         note='descriptive; no verdict; predicted sign = majority sign of the 4 discovery KI-minus-WT differences; labels from the pinned animal manifest'), indent=2) + '\n')


def replay(out_data, cfg=CFG):
    """Independent re-derivation in plain Python (csv module only: no pandas, no numpy, no project functions) of every output row and of the summary. Labels/flags exactly equal, the 12 W difference to 1e-12."""
    def rd(path):
        with open(path, newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
    t3 = rd(cfg['t3']); man = {r['animal_id']: r['genotype'] for r in rd(cfg['animal_manifest'])}
    wt = [a for a in SIX if man[a] == 'WT']; ki = [a for a in SIX if man[a] == 'KI']; saved = rd(Path(out_data) / 'HELDOUT_BY_PC.tsv')
    means = {pop: {r['animal_id']: r for r in rd(Path(cfg['pcmeans']) / f'{pop}_PC_means_by_mouse_12w.tsv')} for pop in TARGETS}
    assert len(wt) == 3 and len(ki) == 3 and len(saved) == sum(len(v) for v in TARGETS.values()); bad = []; seen = [(r['population'], int(r['PC'])) for r in saved]
    assert seen == [(p, k) for p, ks in TARGETS.items() for k in ks]
    for r in saved:
        pop, k = r['population'], int(r['PC']); d = [float(x['diff_KI_minus_WT']) for x in t3 if x['population'] == pop and x['measure'] == f'discovery_PC{k}_mean']
        pos = sum(v > 0 for v in d); neg = sum(v < 0 for v in d); pred = '+' if pos > neg else ('-' if neg > pos else 'none')
        w = [float(means[pop][a][f'PC{k}']) for a in wt]; q = [float(means[pop][a][f'PC{k}']) for a in ki]; diff = sum(q) / 3 - sum(w) / 3
        sign = '+' if diff > 0 else ('-' if diff < 0 else '0'); sep = 'KI_higher' if min(q) > max(w) else ('KI_lower' if max(q) < min(w) else 'none')
        match = 'NA' if pred == 'none' else ('yes' if sign == pred else 'no'); comp = 'NA' if pred == 'none' else ('yes' if {'KI_higher': '+', 'KI_lower': '-', 'none': 'x'}[sep] == pred else 'no')
        ok = (len(d) == 4 and r['predicted_sign'] == pred and r['excluded'] == str(pred == 'none') and abs(float(r['w12_diff_KI_minus_WT']) - diff) < 1e-12 and r['w12_sign'] == sign
              and r['observed_complete_separation'] == sep and r['sign_match'] == match and r['complete_separation_in_predicted_direction'] == comp)
        if not ok: bad.append((pop, k))
    summ = json.loads((Path(out_data) / 'SUMMARY.json').read_text()); sums_ok = True
    for name, sub in [(p, [r for r in saved if r['population'] == p]) for p in TARGETS] + [('both', saved)]:
        pr = [r for r in sub if r['predicted_sign'] != 'none']
        sums_ok = sums_ok and summ[name]['n_with_prediction'] == len(pr) and summ[name]['sign_matches'] == sum(r['sign_match'] == 'yes' for r in pr) and summ[name]['complete_separations_in_predicted_direction'] == sum(r['complete_separation_in_predicted_direction'] == 'yes' for r in pr)
    res = dict(n_rows_checked=len(saved), rows_with_differences=[f'{p}:PC{k}' for p, k in bad], summary_counts_equal=bool(sums_ok)); res['all_pass'] = bool(not bad and sums_ok)
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
    """--check: structure only (no 12 W number is formed): every pre-declared PC has its 4 discovery rows once; the PC-means tables have the PC columns, the six mice and 3 v 3 labels."""
    t3 = pd.read_csv(cfg['t3'], sep='\t'); manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str})
    means_by_pop = {pop: pd.read_csv(Path(cfg['pcmeans']) / f'{pop}_PC_means_by_mouse_12w.tsv', sep='\t', dtype={'animal_id': str}) for pop in TARGETS}
    groups(manifest, means_by_pop)
    for pop, pcs in TARGETS.items():
        for k in pcs: discovery_diffs(t3, pop, k); assert f'PC{k}' in means_by_pop[pop].columns
    return True


def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--check', action='store_true'); g.add_argument('--full', action='store_true'); g.add_argument('--replay', type=Path)
    a = ap.parse_args()
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'): sys.exit('refusing: the checks are asserts; run without -O')
    pins = check_pins(PINS, ROOT)
    if a.check: print(json.dumps(dict(pins_ok=True, n_pins=len(pins), structure_ok=structure_check(CFG)))); return
    if a.replay:
        res = replay(a.replay); docs = ROOT / 'docs' / a.replay.name
        with open(docs / 'REPLAY.json', 'x') as f: f.write(json.dumps(res, indent=2) + '\n')
        print(json.dumps(res)); sys.exit(0 if res['all_pass'] else 1)
    out_data = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_12w_heldout'); out_docs = ROOT / _release_paths.expand('${DATA_ROOT}/phase4_12w_heldout')
    assert not out_data.exists() and not out_docs.exists(), 'outputs exist; never overwrite'
    t0 = time.time(); table, summary, below20 = build(CFG)
    def body():
        write_outputs(out_data, out_docs, table, summary, below20, pins, time.time() - t0); print(json.dumps(dict(summary=summary, mice_below_20_nuclei=below20)))
    guarded(out_docs, body)


if __name__ == '__main__':
    main()

#%% debug
