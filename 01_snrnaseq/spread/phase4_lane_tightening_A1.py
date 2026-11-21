
# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# Phase 4 A1, lane tightening (pre-specified). CODE ONLY until review clears it by hash. DESCRIPTIVE: no verdict, no p-value, no biological reading.
# For the 12 pre-declared genotype PCs (Purkinje 5, 6, 11, 30; Bergmann 2, 4, 5, 6, 8, 9, 20, 29):
#   (1) per mouse: within-mouse SD (ddof = 1) of the nucleus PC scores at 5 / 12 / 18 / 24 / 30 W. Discovery nuclei: the pinned S07 nucleus_PCs.npy (24 mice, IN-SAMPLE: they fitted the discovery PCA);
#       12 W nuclei: the pinned PCs_12w.npy of the accepted full run (6 mice, OUT-OF-SAMPLE projection onto the frozen discovery PCA). Every row carries that label.
#   (2) per age: KI vs WT on the 3 v 3 per-mouse SDs: difference = mean(KI SD) - mean(WT SD); complete separation (strict) KI_higher / KI_lower / none.
#   (3) per genotype: SD vs age trend = the genotype-mean SD at the 5 ages, OLS slope per week on the ages, Spearman rho of mean SD with age (n = 5 ages).
# Labels (WT/KI) and ages come from the pinned animal manifest; the stage A unit table must agree for the 24 discovery mice. A mouse below 20 nuclei keeps its row and carries below_20_nuclei_flag = True.
# Run with: python -I -B ${REPO_ROOT}/01_snrnaseq/spread/phase4_lane_tightening_A1.py --check | --full | --replay DIR
import argparse, csv, hashlib, json, os, statistics, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
SIX = ('2061', '2062', '22018', '2063', '2065', '22019')            # the six 12 W mice
TARGETS = {'Purkinje_cell': (5, 6, 11, 30), 'Bergmann_glia': (2, 4, 5, 6, 8, 9, 20, 29)}
CELLTYPE = {'Purkinje_cell': 'Purkinje cell', 'Bergmann_glia': 'Bergmann glia'}
AGES = (5, 12, 18, 24, 30); FLOOR = 20
IN_SAMPLE = 'in-sample (fitted the discovery PCA)'; OUT_OF_SAMPLE = 'OUT-OF-SAMPLE (12 W, projected onto the frozen discovery PCA)'
CFG = dict(units=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'), core=ROOT / _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01'), r12=ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full'),
           animal_manifest=ROOT / _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'))
PINS = {  # relative path -> full sha256
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'): 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'): 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_PCs.npy'): 'ed535ac9cdb36c51498155641172de92d9b6daa65e611ffa81a705b1a573ea13',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'): '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721',
    _release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Bergmann_glia/nucleus_metadata.tsv.gz'): '23c363d24d24a56f30a8b269f4d3f224c1d9be31ff54ad5f736103e72d0ee540',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Purkinje_cell_PCs_12w.npy'): '3c2eb711b5f3c6777c2f7deb43789337633132194bb09d168544d312ddcfb310',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Bergmann_glia_PCs_12w.npy'): 'c8451ea6a8a4d60607a24a2937a4c47c2c0f2193ce5c35577d1f5c8bfd6fa0d1',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Purkinje_cell_projection_12w.tsv.gz'): '26161c549f37faa58567eca24064603aedec955c5d875996754506e1359d15e4',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/Bergmann_glia_projection_12w.tsv.gz'): '44d6c7eaf54aa361769ff4c7137417526896a93d18f68430c7ff10e26427985c',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/cell_metadata_12w.tsv.gz'): 'f2d3bf5a6987bcd27a2006791e1e2f02d63b6fc9cbb79979b3ac66c4807200fd',
    _release_paths.expand('${DATA_ROOT}/results/phase4_12w_full/SHA256SUMS.txt'): '187e91a7ddaac70c026ee90481cbc02cad6a672be1bed843a01a04c5f9b423fd',
    _release_paths.expand('${DATA_ROOT}/phase4_12w_full/RECORD.json'): '6b0e3c442b1e6ed41ecd2d67ee1d49e59149c22cfcafedadaab7956f0b784123',
    _release_paths.expand('${REPO_ROOT}/data/snrnaseq/sample_table.tsv'): '7c2f2d4d41ae6e06b20f8c18071110ae776d93b0b46978f510e19044b912bb07'}


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


def load_population(cfg, pop, manifest):
    """Nucleus PC scores of one population: the 24 discovery mice (in-sample) then the six 12 W mice (out-of-sample); the mouse of every row is found by obs_name / unit_index, never by position alone."""
    name = f'discovery__{pop}'; nm = pd.read_csv(Path(cfg['core']) / name / 'nucleus_metadata.tsv.gz', sep='\t'); pcs = np.load(Path(cfg['core']) / name / 'nucleus_PCs.npy')
    units = pd.read_csv(cfg['units'], sep='\t', dtype={'animal_id': str}); assert (units.unit_index.to_numpy() == np.arange(len(units))).all()
    u = units.iloc[nm.unit_index.to_numpy()]; g = manifest.set_index('animal_id')
    assert len(nm) == len(pcs) and (u.cell_type_raw == CELLTYPE[pop]).all(), pop
    assert (u.genotype.to_numpy() == g.genotype.loc[u.animal_id].to_numpy()).all() and (u.age_weeks.to_numpy() == g.age_weeks.loc[u.animal_id].to_numpy()).all(), ('stage A unit table disagrees with the animal manifest', pop)
    p12 = np.load(Path(cfg['r12']) / f'{pop}_PCs_12w.npy'); proj = pd.read_csv(Path(cfg['r12']) / f'{pop}_projection_12w.tsv.gz', sep='\t', dtype={'animal_id': str}, usecols=['obs_name', 'animal_id', 'age_weeks'])
    cells = pd.read_csv(Path(cfg['r12']) / 'cell_metadata_12w.tsv.gz', sep='\t', dtype={'animal_id': str}); c = cells[cells.cell_type_raw == CELLTYPE[pop]]
    assert len(proj) == len(p12) and proj.obs_name.tolist() == c.obs_name.tolist() and proj.animal_id.tolist() == c.animal_id.tolist(), ('12 W PC rows do not map to the cell table', pop)
    assert set(proj.animal_id) == set(SIX) and (proj.age_weeks == 12).all() and (g.age_weeks.loc[list(SIX)] == 12).all()
    ids = np.concatenate([u.animal_id.to_numpy().astype(object), proj.animal_id.to_numpy().astype(object)]); roles = np.array([IN_SAMPLE] * len(u) + [OUT_OF_SAMPLE] * len(proj), dtype=object)
    return np.vstack([pcs, p12]), ids, roles


def sd_rows(pop, pcs, ids, roles, manifest):
    """Per mouse x PC: n nuclei and the within-mouse SD (ddof = 1) of the nucleus scores of PC k (column k - 1)."""
    g = manifest.set_index('animal_id'); rows = []
    for k in TARGETS[pop]:
        col = pcs[:, k - 1]
        for a in sorted(set(ids.tolist())):
            m = ids == a; assert len(set(roles[m].tolist())) == 1
            rows.append(dict(population=pop, PC=k, animal_id=a, genotype=g.genotype[a], age_weeks=int(g.age_weeks[a]), n_nuclei=int(m.sum()), below_20_nuclei_flag=bool(m.sum() < FLOOR), sd=float(np.std(col[m], ddof=1)), role=roles[m][0]))
    return rows


def by_age(sd):
    rows = []
    for (pop, k, age), grp in sd.groupby(['population', 'PC', 'age_weeks'], sort=True):
        w = grp[grp.genotype == 'WT'].sd.to_numpy(); q = grp[grp.genotype == 'KI'].sd.to_numpy(); assert len(w) == 3 and len(q) == 3, (pop, k, age)
        sep = 'KI_higher' if q.min() > w.max() else ('KI_lower' if q.max() < w.min() else 'none')
        rows.append(dict(population=pop, PC=int(k), age_weeks=int(age), role=grp.role.iloc[0], n_WT=3, n_KI=3, WT_mean_sd=float(w.mean()), KI_mean_sd=float(q.mean()), diff_KI_minus_WT=float(q.mean() - w.mean()), separation=sep))
    return pd.DataFrame(rows)


def trend(sd):
    rows = []; x = np.array(AGES, float)
    for (pop, k, geno), grp in sd.groupby(['population', 'PC', 'genotype'], sort=True):
        m = grp.groupby('age_weeks').sd.mean(); assert list(m.index) == list(AGES), (pop, k, geno)
        y = m.to_numpy(); slope = float(((x - x.mean()) * (y - y.mean())).sum() / ((x - x.mean()) ** 2).sum()); rho = float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])
        row = dict(population=pop, PC=int(k), genotype=geno); row.update({f'mean_sd_{a}W': float(v) for a, v in zip(AGES, y)}); row.update(slope_sd_per_week=slope, spearman_rho_sd_vs_age=rho); rows.append(row)
    return pd.DataFrame(rows)


def summarize(ba, tr):
    out = {}
    for name, b, t in [(p, ba[ba.population == p], tr[tr.population == p]) for p in TARGETS] + [('both', ba, tr)]:
        per_age = {}
        for a in AGES:
            s = b[b.age_weeks == a]; n = int(len(s))
            per_age[f'{a}W'] = dict(n_tests=n, diff_KI_higher=int((s.diff_KI_minus_WT > 0).sum()), diff_KI_lower=int((s.diff_KI_minus_WT < 0).sum()), complete_separation_KI_higher=int((s.separation == 'KI_higher').sum()),
                                    complete_separation_KI_lower=int((s.separation == 'KI_lower').sum()), chance_expected_separations_each_direction=n * 0.05)
        out[name] = dict(by_age=per_age, trend={geno: dict(n_pcs=int(len(t[t.genotype == geno])), n_negative_slope=int((t[t.genotype == geno].slope_sd_per_week < 0).sum()), n_positive_slope=int((t[t.genotype == geno].slope_sd_per_week > 0).sum())) for geno in ('WT', 'KI')})
    return out


def build(cfg):
    manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}); rows = []
    for pop in TARGETS:
        pcs, ids, roles = load_population(cfg, pop, manifest); rows += sd_rows(pop, pcs, ids, roles, manifest)
    sd = pd.DataFrame(rows); ba = by_age(sd); tr = trend(sd)
    return sd, ba, tr, summarize(ba, tr)


def write_outputs(out_data, out_docs, sd, ba, tr, summary, pins, seconds):
    out_data.mkdir(parents=True); out_docs.mkdir(parents=True)
    sd.to_csv(out_data / 'SD_BY_MOUSE_AGE_PC.tsv', sep='\t', index=False); ba.to_csv(out_data / 'KI_VS_WT_BY_AGE_PC.tsv', sep='\t', index=False); tr.to_csv(out_data / 'TREND_BY_GENOTYPE_PC.tsv', sep='\t', index=False)
    (out_data / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n'); (out_data / 'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(out_data.iterdir()) if p.is_file()))
    flagged = sd[sd.below_20_nuclei_flag].drop_duplicates(['population', 'animal_id'])[['population', 'animal_id', 'n_nuclei']].to_dict('records')
    (out_docs / 'RECORD.json').write_text(json.dumps(dict(seconds=seconds, summary=summary, mice_below_20_nuclei=flagged, pins=pins, script_sha256=sha(Path(__file__)),
                                                         note='descriptive; no verdict; discovery rows are in-sample, 12 W rows are out-of-sample; SD uses ddof = 1; labels and ages from the pinned animal manifest'), indent=2) + '\n')


def replay(out_data, cfg=CFG):
    """Independent re-derivation: every per-mouse SD again with statistics.stdev on plain Python floats (numpy only loads the arrays); the by-age and trend tables again from the SAVED SD table in plain Python (csv module, manual average ranks).
    Labels/flags exactly equal, floats to 1e-9 (relative)."""
    def rd(path):
        with open(path, newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
    def close(a, b): return abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b)))
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]: j += 1
            for t in range(i, j + 1): r[order[t]] = (i + j) / 2 + 1
            i = j + 1
        return r
    out_data = Path(out_data); man = {r['animal_id']: r for r in rd(cfg['animal_manifest'])}; saved = rd(out_data / 'SD_BY_MOUSE_AGE_PC.tsv'); bad = []
    expected_n = sum(len(v) for v in TARGETS.values()) * 30; assert len(saved) == expected_n, len(saved)
    cache = {}
    for pop in TARGETS:
        name = f'discovery__{pop}'; nm = pd.read_csv(Path(cfg['core']) / name / 'nucleus_metadata.tsv.gz', sep='\t')
        units = rd(cfg['units']); ids = [units[int(i)]['animal_id'] for i in nm.unit_index]; pcs = np.load(Path(cfg['core']) / name / 'nucleus_PCs.npy'); p12 = np.load(Path(cfg['r12']) / f'{pop}_PCs_12w.npy')
        proj = rd_gz(Path(cfg['r12']) / f'{pop}_projection_12w.tsv.gz'); cache[pop] = (ids + [r['animal_id'] for r in proj], np.vstack([pcs, p12]))
    for r in saved:
        pop, k, a = r['population'], int(r['PC']), r['animal_id']; ids, pcs = cache[pop]; vals = [float(pcs[i, k - 1]) for i in range(len(ids)) if ids[i] == a]
        ok = (len(vals) == int(r['n_nuclei']) and close(statistics.stdev(vals), r['sd']) and r['genotype'] == man[a]['genotype'] and int(r['age_weeks']) == int(float(man[a]['age_weeks'])) and r['below_20_nuclei_flag'] == str(len(vals) < FLOOR)
              and r['role'] == (OUT_OF_SAMPLE if a in SIX else IN_SAMPLE))
        if not ok: bad.append((pop, k, a))
    ba_saved = rd(out_data / 'KI_VS_WT_BY_AGE_PC.tsv'); bad_age = []
    for r in ba_saved:
        pop, k, age = r['population'], int(r['PC']), int(r['age_weeks'])
        w = [float(s['sd']) for s in saved if s['population'] == pop and int(s['PC']) == k and int(s['age_weeks']) == age and s['genotype'] == 'WT']
        q = [float(s['sd']) for s in saved if s['population'] == pop and int(s['PC']) == k and int(s['age_weeks']) == age and s['genotype'] == 'KI']
        sep = 'KI_higher' if min(q) > max(w) else ('KI_lower' if max(q) < min(w) else 'none')
        if not (len(w) == 3 and len(q) == 3 and close(r['diff_KI_minus_WT'], sum(q) / 3 - sum(w) / 3) and r['separation'] == sep and close(r['WT_mean_sd'], sum(w) / 3) and close(r['KI_mean_sd'], sum(q) / 3)): bad_age.append((pop, k, age))
    tr_saved = rd(out_data / 'TREND_BY_GENOTYPE_PC.tsv'); bad_tr = []
    for r in tr_saved:
        pop, k, geno = r['population'], int(r['PC']), r['genotype']
        y = [statistics.fmean(float(s['sd']) for s in saved if s['population'] == pop and int(s['PC']) == k and s['genotype'] == geno and int(s['age_weeks']) == a) for a in AGES]; x = [float(a) for a in AGES]
        mx, my = sum(x) / 5, sum(y) / 5; slope = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y)) / sum((xi - mx) ** 2 for xi in x)
        rx, ry = ranks(x), ranks(y); mrx, mry = sum(rx) / 5, sum(ry) / 5; den = (sum((a - mrx) ** 2 for a in rx) * sum((b - mry) ** 2 for b in ry)) ** 0.5
        rho = sum((a - mrx) * (b - mry) for a, b in zip(rx, ry)) / den if den else float('nan')
        if not (close(r['slope_sd_per_week'], slope) and (close(r['spearman_rho_sd_vs_age'], rho) if den else r['spearman_rho_sd_vs_age'] in ('', 'nan')) and all(close(r[f'mean_sd_{a}W'], v) for a, v in zip(AGES, y))): bad_tr.append((pop, k, geno))
    summ = json.loads((out_data / 'SUMMARY.json').read_text()); sums_ok = True
    for name, pops in [(p, [p]) for p in TARGETS] + [('both', list(TARGETS))]:
        for a in AGES:
            sub = [r for r in ba_saved if r['population'] in pops and int(r['age_weeks']) == a]; e = summ[name]['by_age'][f'{a}W']
            sums_ok = sums_ok and e['n_tests'] == len(sub) and e['diff_KI_higher'] == sum(float(r['diff_KI_minus_WT']) > 0 for r in sub) and e['diff_KI_lower'] == sum(float(r['diff_KI_minus_WT']) < 0 for r in sub) \
                and e['complete_separation_KI_higher'] == sum(r['separation'] == 'KI_higher' for r in sub) and e['complete_separation_KI_lower'] == sum(r['separation'] == 'KI_lower' for r in sub)
    res = dict(n_sd_rows_checked=len(saved), sd_rows_with_differences=[list(b) for b in bad], by_age_rows_with_differences=[list(b) for b in bad_age], trend_rows_with_differences=[list(b) for b in bad_tr], summary_counts_equal=bool(sums_ok)); res['all_pass'] = bool(not bad and not bad_age and not bad_tr and sums_ok)
    return res


def rd_gz(path):
    import gzip
    with gzip.open(path, 'rt', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))


def guarded(out_docs, fn):
    """Run fn; on ANY failure after the output dirs exist, exclusive-create FAILED.json (traceback) in out_docs, then re-raise (no retry, nothing deleted)."""
    try: return fn()
    except BaseException:
        import traceback
        if Path(out_docs).exists():
            with open(Path(out_docs) / 'FAILED.json', 'x') as f: f.write(json.dumps(dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), traceback=traceback.format_exc()), indent=2) + '\n')
        raise


def structure_check(cfg):
    """--check: structure only (no SD is formed): every population loads, the row mapping asserts hold, every pre-declared PC column exists, each of the 5 ages has 3 WT and 3 KI mice."""
    manifest = pd.read_csv(cfg['animal_manifest'], sep='\t', dtype={'animal_id': str}); g = manifest.set_index('animal_id')
    for pop in TARGETS:
        pcs, ids, roles = load_population(cfg, pop, manifest); assert pcs.shape[1] >= max(TARGETS[pop]); mice = sorted(set(ids.tolist()))
        assert len(mice) == 30 and sorted((int(g.age_weeks[a]), g.genotype[a]) for a in mice) == sorted((age, geno) for age in AGES for geno in ('WT', 'KI') for _ in range(3)), pop
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
    out_data = ROOT / _release_paths.expand('${DATA_ROOT}/results/phase4_lane_tightening_A1'); out_docs = ROOT / _release_paths.expand('${DATA_ROOT}/phase4_lane_tightening_A1')
    assert not out_data.exists() and not out_docs.exists(), 'outputs exist; never overwrite'
    t0 = time.time(); sd, ba, tr, summary = build(CFG)
    def body():
        write_outputs(out_data, out_docs, sd, ba, tr, summary, pins, time.time() - t0); print(json.dumps(summary))
    guarded(out_docs, body)


if __name__ == '__main__':
    main()

#%% debug
