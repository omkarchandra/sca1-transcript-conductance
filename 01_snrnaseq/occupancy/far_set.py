"""Far-set KI/WT ratios, Purkinje cells .  EXPLORATORY; KI nuclei selected on the genotype lanes; ratios biased upward by construction.

Stage A (python): lane direction per age, far-set selection, SELECTION.tsv, OCCUPANCY.tsv, OCCUPANCY_SUMMARY.md, pseudobulk columns.
Stage B (R, far_set_fit.R): control (all KI nuclei) first, gate against the S08 release, then the two far-set rules.
python -I -B 01_snrnaseq/occupancy/far_set.py --mode real --out DIR --rscript R --r-script F.R --r-sha256 H --cleared-sha256 H
       python -I -B 01_snrnaseq/occupancy/far_set.py --mode fixture --inputs DIR --out DIR --rscript R --r-script F.R --r-sha256 H
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

import argparse, gzip, hashlib, io, json, math, os, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
DISCLAIMER = 'Exploratory; KI nuclei selected on the genotype lanes; ratios biased upward by construction.'
AGES_ALL = (5, 18, 24, 30)
AGES = (18, 24, 30)
LANE_CAND = (5, 6, 11, 30)
EXPECTED_LANES = {5: (5, 6, 11), 18: (5, 6, 11, 30), 24: (5, 6, 11, 30), 30: (5, 6, 30)}   # from pc_check_by_pc.tsv (separating PCs); re-derived and asserted at run time
MODEL_GENES = ['Kcnc3', 'Kcnma1', 'Kcnd3', 'Cacna1a', 'Cacna1g', 'Slc8a1', 'Itpr1', 'Atp2a3']
CUTS = (0.50, 0.90, 0.95)
FLOOR = 20
GATE_EIGHT_MAX = 0.02      # hard stop: control vs S08 release per-age log2FC, eight genes
GATE_ALL_MAX = 0.10        # hard stop: control vs S08 release per-age log2FC, any core gene
T3_DIGITS = 6              # the pinned T3 table stores per-mouse means rounded to 6 significant digits; hard stop if a mean differs by more than half a unit of the last digit
REAL_PINS = {
    'nucleus_pcs':  (_release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_PCs.npy'), 'f2a6ed3a7f90068228f675eb8e17f8d41860402d8d8c2df12abb8cabbfe5ec3b'),
    'nucleus_meta': (_release_paths.expand('${DATA_ROOT}/results/s07_core_representation_v01/discovery__Purkinje_cell/nucleus_metadata.tsv.gz'), '0f02b61d45667f203d410dce7751bf80c7c78a00a6af65afc30ba97039317721'),
    'unit_meta':    (_release_paths.expand('${REPO_ROOT}/data/snrnaseq/unit_metadata.tsv'), 'd2a0b8b202d8eb5823ec6c8f1907630a0f441352dbe0b4ed63ca9d918a47d657'),
    't3_table':     (_release_paths.expand('${REPO_ROOT}/results/s12_lanes/kivswt_within_age.tsv'), '01f29ce2b9c9e648b37e478ed408e9a3d276df5240d776e6996fa1cd24565914'),
    'pc_check':     (_release_paths.expand('${REPO_ROOT}/results/s12_lanes/pc_check_by_pc.tsv'), '27701bd432255403c50f399ad9a9f98a374b8cdccdfb2e8c22e87849b3bc4b31'),
    'a3_lane':      (_release_paths.expand('${DATA_ROOT}/results/phase4_a3_bimodality/PER_MOUSE_LANE.tsv'), '48931e8971ba8b137dce6532c2af34fe18a9ae98c453166b2a06947ee1fa43e6'),
    's08_counts':   (_release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/raw_counts.tsv.gz'), '3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778'),
    's08_unit_meta': (_release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/unit_metadata.tsv'), '11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4'),
    's08_gene_meta': (_release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/gene_metadata.tsv'), '765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3'),
    'release':      (_release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz'), '40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298'),
    'expr_rules':   (_release_paths.expand('${REPO_ROOT}/01_snrnaseq/s08_limma/s08_expression_rules.R'), 'e1b64c9f6259812cc9d758cea6867f3ef8f4150e1bc2553a46e33d5660449981'),
    's08_limma_ref': (_release_paths.expand('${REPO_ROOT}/01_snrnaseq/s08_limma/s08_limma_real.R'), 'e77fea81213ed9dec2dcdf97503e6999f0ca555b4000167f9dde8f830cf0143c'),
    'h5ad':         (_release_paths.expand('${DATA_ROOT}/data/adata_leiden_0.175_annotated_cell_types_v3.h5ad'), '9fb8013c5879588d04008a8778a15a31708b748bd7dad74c9ea85fae140dad86'),
}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def wtext(path, text):
    with open(path, 'x') as f:
        f.write(text)


def wtsv(df, path, comment=True, float_format=None):
    buf = io.StringIO()
    if comment:
        buf.write('# ' + DISCLAIMER + '\n')
    df.to_csv(buf, sep='\t', index=False, float_format=float_format)
    if str(path).endswith('.gz'):
        with gzip.open(path, 'xt') as f:
            f.write(buf.getvalue())
    else:
        wtext(path, buf.getvalue())


def resolve_pins(args):
    if args.mode == 'real':
        pins = {k: (ROOT / rel, h) for k, (rel, h) in REAL_PINS.items()}
    else:
        base = Path(args.inputs)
        pj = json.loads((base / 'PINS.json').read_text())
        pins = {k: (base / rel, h) for k, (rel, h) in pj.items()}
        assert set(pins) == set(REAL_PINS), 'fixture pins must name the same inputs as the real pins'
    out = {}
    for k, (p, h) in pins.items():
        assert p.exists(), f'pin missing: {k} {p}'
        if k == 'h5ad':
            out[k] = (p, h)          # hashed separately (before and after the read)
            continue
        assert sha(p) == h, f'sha256 mismatch for pin: {k}'
        out[k] = (p, h)
    return out


def parse_values(s):
    d = {}
    for part in str(s).split(';'):
        k, v = part.split(':')
        d[k.strip()] = float(v)
    return d


def stage_a(args, pins, out, rec):
    P = {k: v[0] for k, v in pins.items()}
    X = np.load(P['nucleus_pcs'])
    nm = pd.read_csv(P['nucleus_meta'], sep='\t')
    um = pd.read_csv(P['unit_meta'], sep='\t', dtype={'animal_id': str})
    assert X.ndim == 2 and X.shape[0] == len(nm) and X.shape[1] >= 30 and np.isfinite(X).all()
    assert nm.obs_name.is_unique and (um.index.to_numpy() == um.unit_index.to_numpy()).all()
    pk = um[um.cell_type_raw == 'Purkinje cell']
    assert len(pk) == 24 and set(pk.age_weeks) == set(AGES_ALL) and (pk.groupby(['age_weeks', 'genotype']).size() == 3).all()
    assert set(nm.unit_index) == set(pk.unit_index) and (nm.groupby('unit_index').size().reindex(pk.unit_index).to_numpy() == pk.n_nuclei.to_numpy()).all()
    cell = nm.merge(pk[['unit_index', 'animal_id', 'age_weeks', 'genotype', 'sex', 'n_nuclei']], on='unit_index', how='left', validate='many_to_one')
    assert cell.animal_id.notna().all() and (cell.obs_name.to_numpy() == nm.obs_name.to_numpy()).all()
    cell['row'] = np.arange(len(cell))
    mice = pk[['animal_id', 'age_weeks', 'genotype', 'sex', 'n_nuclei', 'unit_index']].reset_index(drop=True)

    # lane PCs per age from the pinned PC check table
    pcc = pd.read_csv(P['pc_check'], sep='\t')
    pcc = pcc[pcc.population == 'Purkinje_cell'].set_index('PC')
    lanes, sign = {}, {}
    for a in AGES_ALL:
        col = f'sep_{a}W'
        lanes[a] = tuple(k for k in LANE_CAND if pcc.loc[k, col] != 'none')
        sign[a] = {k: (1 if pcc.loc[k, col] == 'KI_higher' else -1) for k in lanes[a]}
        assert lanes[a] == EXPECTED_LANES[a], ('lane PCs differ from the pre-registered set', a, lanes[a])
    # per-mouse mean PC scores, cross-checked against the pinned T3 table
    t3 = pd.read_csv(P['t3_table'], sep='\t', dtype={'WT_values': str, 'KI_values': str})
    t3 = t3[t3.population == 'Purkinje_cell']
    mean_diff, maxdev = {}, 0.0
    for a in AGES_ALL:
        am = mice[mice.age_weeks == a]
        ki = am[am.genotype == 'KI'].animal_id.tolist(); wt = am[am.genotype == 'WT'].animal_id.tolist()
        d = []
        for k in LANE_CAND:
            mm = {m: X[(cell.animal_id == m).to_numpy(), k - 1].mean() for m in ki + wt}
            row = t3[(t3.age_weeks == a) & (t3.measure == f'discovery_PC{k}_mean')]
            assert len(row) == 1, ('T3 row', a, k)
            row = row.iloc[0]
            ref = {**parse_values(row.WT_values), **parse_values(row.KI_values)}
            assert set(ref) == set(mm), ('T3 mouse ids', a, k)
            for m in mm:
                e = math.floor(math.log10(abs(ref[m]))) if ref[m] != 0 else 0
                tol = 0.5 * 10.0 ** (e - (T3_DIGITS - 1)) * 1.001
                dev = abs(ref[m] - mm[m]); maxdev = max(maxdev, dev)
                assert dev <= tol, ('per-mouse PC mean differs from the pinned T3 table', a, k, m, dev, tol)
            diff = np.mean([mm[m] for m in ki]) - np.mean([mm[m] for m in wt])
            if k in lanes[a]:
                assert np.sign(diff) == sign[a][k], ('lane sign disagrees with pc_check', a, k)
                lo, hi = sorted([mm[m] for m in ki]), sorted([mm[m] for m in wt])
                assert (lo[0] > hi[-1]) if sign[a][k] == 1 else (lo[-1] < hi[0]), ('lane PC is not completely separated', a, k)
            d.append(diff)
        mean_diff[a] = np.array(d)
    rec['t3_max_deviation_of_per_mouse_means'] = maxdev

    # lane direction and scores
    score = np.full(len(cell), np.nan); direction = {}
    for a in AGES_ALL:
        idx = [LANE_CAND.index(k) for k in lanes[a]]
        dvec = mean_diff[a][idx]
        u = dvec / np.linalg.norm(dvec)
        direction[a] = {f'PC{k}': float(v) for k, v in zip(lanes[a], u)}
        sel = (cell.age_weeks == a).to_numpy()
        score[sel] = X[sel][:, [k - 1 for k in lanes[a]]] @ u
    assert np.isfinite(score).all()
    cell['score'] = score
    rec['lane_pcs'] = {str(a): list(lanes[a]) for a in AGES_ALL}
    rec['lane_unit_vector'] = {str(a): direction[a] for a in AGES_ALL}
    for a in AGES_ALL:
        s = cell[cell.age_weeks == a]
        assert s[s.genotype == 'KI'].groupby('animal_id').score.mean().mean() > s[s.genotype == 'WT'].groupby('animal_id').score.mean().mean()

    # thresholds and selection
    thr = {}
    for a in AGES_ALL:
        w = cell[(cell.age_weeks == a) & (cell.genotype == 'WT')].score.to_numpy()
        thr[a] = {c: float(np.quantile(w, c)) for c in CUTS}
    rec['wt_quantile_thresholds'] = {str(a): {str(c): v for c, v in thr[a].items()} for a in AGES_ALL}
    cell['sel_wt90'] = False; cell['sel_top_third'] = False
    third_cut = {}
    for a in AGES:
        ki = cell[(cell.age_weeks == a) & (cell.genotype == 'KI')]
        cell.loc[ki.index[ki.score.to_numpy() > thr[a][0.90]], 'sel_wt90'] = True
        order = np.argsort(-ki.score.to_numpy(), kind='stable')
        kk = int(math.ceil(len(ki) / 3))
        cell.loc[ki.index[order[:kk]], 'sel_top_third'] = True
        third_cut[a] = float(ki.score.to_numpy()[order[kk - 1]])
    sel_rows = []
    for rule, col in (('wt90', 'sel_wt90'), ('top_third', 'sel_top_third')):
        for a in AGES:
            for _, m in mice[(mice.age_weeks == a) & (mice.genotype == 'KI')].iterrows():
                c = cell[cell.animal_id == m.animal_id]
                n = int(c[col].sum())
                sel_rows.append(dict(age=a, mouse=m.animal_id, sex=m.sex, rule=rule, n_nuclei=len(c), n_selected=n, share=n / len(c),
                                     threshold=thr[a][0.90] if rule == 'wt90' else third_cut[a], below_s08_floor_20=n < FLOOR))
    sel = pd.DataFrame(sel_rows)
    assert (sel.n_selected >= 1).all(), ('a KI mouse has no far-set nucleus under a rule: STOP and report', sel[sel.n_selected < 1].to_dict('records'))
    wtsv(sel, out / 'SELECTION.tsv', float_format='%.6g')

    # OCCUPANCY
    a3 = pd.read_csv(P['a3_lane'], sep='\t', dtype={'animal_id': str})
    a3 = a3[a3.population.astype(str).str.startswith('Purkinje')]
    assert len(a3) == 48, len(a3)
    occ = []
    for _, m in mice.iterrows():
        c = cell[cell.animal_id == m.animal_id]; a = int(m.age_weeks); s = c.score.to_numpy()
        r = dict(age=a, genotype=m.genotype, mouse=m.animal_id, n_nuclei_total=len(c), lane_score_median=float(np.median(s)),
                 lane_score_q25=float(np.quantile(s, .25)), lane_score_q75=float(np.quantile(s, .75)), lane_score_iqr=float(np.quantile(s, .75) - np.quantile(s, .25)))
        for q in CUTS:
            n = int((s > thr[a][q]).sum()); r[f'n_beyond_WT_p{round(q*100)}'] = n; r[f'share_beyond_WT_p{round(q*100)}'] = n / len(s)
        for k in LANE_CAND:
            g = a3[(a3.animal_id == m.animal_id) & (a3.PC == k)]
            r[f'a3_call_PC{k}'] = g.two_group_call.iloc[0] if len(g) == 1 else ''
        occ.append(r)
    occ = pd.DataFrame(occ)
    wtsv(occ, out / 'OCCUPANCY.tsv', float_format='%.6g')
    lines = ['# OCCUPANCY_SUMMARY', '', DISCLAIMER, '',
             'Lane score = projection on the unit vector of the KI-WT difference of per-mouse mean scores over the separating Purkinje lane PCs at that age (5 W is occupancy only, no ratios).',
             'Shares are the mean of the 3 mice per group; "beyond" = score strictly greater than that percentile of the pooled WT nuclei of the same age.', '']
    for a in AGES_ALL:
        lines += [f'## {a} weeks (lane PCs {list(lanes[a])})', '', '| group | n mice | mean share beyond WT p50 | p90 | p95 |', '|---|---|---|---|---|']
        for g in ('KI', 'WT'):
            o = occ[(occ.age == a) & (occ.genotype == g)]
            lines.append(f'| {g} | {len(o)} | ' + ' | '.join(f'{o[f"share_beyond_WT_p{round(q*100)}"].mean():.3f}' for q in CUTS) + ' |')
        lines.append('')
        if a in AGES:
            for rule in ('wt90', 'top_third'):
                sr = sel[(sel.age == a) & (sel.rule == rule)]
                lines.append(f'Far-set nuclei per KI mouse, rule {rule}: ' + ', '.join(f'{r.mouse}: {r.n_selected} of {r.n_nuclei}' + (' (below 20)' if r.below_s08_floor_20 else '') for r in sr.itertuples()))
            lines.append('')
    wtext(out / 'OCCUPANCY_SUMMARY.md', '\n'.join(lines) + '\n')
    return cell, mice, sel, thr


def read_counts_rows(h5ad, names_wanted):
    import h5py
    from anndata.io import read_elem
    with h5py.File(h5ad, 'r', rdcc_nbytes=1 << 28) as f:
        var = read_elem(f['var'])
        key = f['obs'].attrs['_index']; key = key.decode() if isinstance(key, bytes) else key
        obj = f['obs'][key]
        names = pd.Index(np.asarray(obj.asstr()[...] if isinstance(obj, h5py.Dataset) else read_elem(obj), dtype=str))
        assert names.is_unique
        pos = names.get_indexer(names_wanted)
        assert (pos >= 0).all(), 'nucleus names missing from the h5ad'
        g = f['layers/counts']; ptr = g['indptr'][...].astype(np.int64)
        import scipy.sparse as sp
        order = np.argsort(pos)
        datas, inds, lens = [None] * len(pos), [None] * len(pos), np.zeros(len(pos), dtype=np.int64)
        for i in order:
            lo, hi = int(ptr[pos[i]]), int(ptr[pos[i] + 1])
            datas[i] = g['data'][lo:hi].astype(np.float64); inds[i] = g['indices'][lo:hi]; lens[i] = hi - lo
        indptr = np.concatenate([[0], np.cumsum(lens)])
        M = sp.csr_matrix((np.concatenate(datas), np.concatenate(inds), indptr), shape=(len(pos), len(var)))
    return M, var


def stage_pseudobulk(pins, out, cell, mice, rec):
    import scipy.sparse as sp
    P = {k: v[0] for k, v in pins.items()}
    h0 = sha(P['h5ad']); assert h0 == pins['h5ad'][1], 'h5ad sha256 mismatch'
    M, var = read_counts_rows(P['h5ad'], pd.Index(cell.obs_name))
    assert sha(P['h5ad']) == h0, 'h5ad changed during the read'
    gene_ids = var.gene_ids.to_numpy().astype(str)
    gm = pd.read_csv(P['s08_gene_meta'], sep='\t')
    assert (gm.stable_id.to_numpy().astype(str) == gene_ids).all(), 'gene order differs from the S08 gene metadata'
    assert (M.data >= 0).all() and (M.data == np.floor(M.data)).all()
    rs = np.asarray(M.sum(axis=1)).ravel()
    assert np.array_equal(rs, cell.full_count_total.to_numpy().astype(float)), 'nucleus row sums differ from full_count_total'
    mouse_index = {m: i for i, m in enumerate(mice.animal_id)}
    mi = cell.animal_id.map(mouse_index).to_numpy()

    def sums(mask):
        ind = sp.csr_matrix((np.ones(int(mask.sum())), (mi[mask], np.flatnonzero(mask))), shape=(len(mice), len(cell)))
        return np.asarray((ind @ M).todense()).astype(np.int64)
    base = sums(np.ones(len(cell), dtype=bool))
    # exact match with the S08 pseudobulk columns of these units
    s08 = pd.read_csv(P['s08_counts'], sep='\t', index_col=0)
    um8 = pd.read_csv(P['s08_unit_meta'], sep='\t', dtype={'animal_id': str})
    assert um8.unit_label.tolist() == [f'U{i:03d}' for i in range(len(um8))], 'S08 unit labels are not U000..'
    cols = [f'U{int(u):03d}' for u in mice.unit_index]
    assert all(um8.animal_id.iloc[int(u)] == m and um8.cell_type_raw.iloc[int(u)] == 'Purkinje cell' for u, m in zip(mice.unit_index, mice.animal_id)), 'S08 unit rows do not match the mice'
    assert (s08.index.to_numpy().astype(str) == gene_ids).all()
    assert np.array_equal(s08[cols].to_numpy().astype(np.int64), base.T), 'pseudobulk of all nuclei differs from the S08 raw_counts columns'
    rec['pseudobulk_matches_s08_raw_counts_exactly'] = True
    variants = {'control': base.copy(), 'wt90': base.copy(), 'top_third': base.copy()}
    ncols = {v: mice.n_nuclei.to_numpy().copy() for v in variants}
    for rule, col in (('wt90', 'sel_wt90'), ('top_third', 'sel_top_third')):
        for i, m in mice.iterrows():
            if m.age_weeks in AGES and m.genotype == 'KI':
                mask = ((cell.animal_id == m.animal_id) & cell[col]).to_numpy()
                variants[rule][i] = sums(mask)[i]
                ncols[rule][i] = int(mask.sum())
                assert (variants[rule][i] <= base[i]).all()
    pb = out / 'pseudobulk'; pb.mkdir()
    meta_rows = []
    for v, mat in variants.items():
        df = pd.DataFrame(mat.T, index=gene_ids, columns=[f'M{m}' for m in mice.animal_id]); df.index.name = 'stable_id'
        with gzip.open(pb / f'pb_{v}.tsv.gz', 'xt') as f:
            df.to_csv(f, sep='\t')
        for i, m in mice.iterrows():
            meta_rows.append(dict(variant=v, column=f'M{m.animal_id}', animal_id=m.animal_id, age_weeks=int(m.age_weeks), genotype=m.genotype, sex=m.sex, n_nuclei_in_column=int(ncols[v][i])))
    pd.DataFrame(meta_rows).to_csv(pb / 'pb_columns.tsv', sep='\t', index=False)
    return gm


def run_r(args, pins, variants, indir, outdir):
    P = {k: v[0] for k, v in pins.items()}
    cmd = [args.rscript, args.r_script, f'--indir={indir}', f'--outdir={outdir}', f'--rules={P["expr_rules"]}', f'--variants={",".join(variants)}', f'--floor={FLOOR}']
    r = subprocess.run(cmd, capture_output=True, text=True)
    wtext(Path(outdir).parent / f'R_console_{Path(outdir).name}.txt', r.stdout + '\n---stderr---\n' + r.stderr)
    assert r.returncode == 0, ('R stage failed', r.stderr[-2000:])


def load_fit(fitdir):
    idx = pd.read_csv(Path(fitdir) / 'FIT_INDEX.tsv', sep='\t', dtype={'ki_mice': str})
    return idx


def table(fitdir, row):
    return pd.read_csv(Path(fitdir) / 'tables' / row.table_file, sep='\t') if isinstance(row.table_file, str) else None


def gate(fitdir, pins, gm, rec):
    P = {k: v[0] for k, v in pins.items()}
    rel = pd.read_csv(P['release'], sep='\t', low_memory=False)
    rel = rel[(rel.population == 'purkinje') & (rel.family == 'core_primary')]
    idx = load_fit(fitdir)
    eight = gm[gm.dataset_symbol.isin(MODEL_GENES)]
    assert len(eight) == 8 and eight.dataset_symbol.is_unique
    rows, e_max, a_max = [], 0.0, 0.0
    for a in AGES:
        r = rel[rel.age_weeks == a].set_index('stable_id')
        for fit, relcol in (('perage', 'fallback_perage_log2FC'), ('pooled', 'log2FC')):
            ix = idx[(idx.variant == 'control') & (idx.floor == 1) & (idx.fit == fit) & (idx.age == a)]
            assert len(ix) == 1 and ix.iloc[0].status == 'FITTED', ('control fit missing', a, fit)
            t = table(fitdir, ix.iloc[0]).set_index('stable_id')
            common = r.index.intersection(t.index)
            assert len(common) >= 0.9 * len(r), ('control table covers too few release genes', a, fit, len(common), len(r))
            d = (t.loc[common, 'logFC'] - r.loc[common, relcol]).abs()
            de = d.reindex(eight.stable_id).dropna()
            assert len(de) == 8, ('eight genes not all in the control table', a, fit)
            e_max = max(e_max, float(de.max())); a_max = max(a_max, float(d.max()))
            rows.append(dict(age=a, fit=fit, n_common=len(common), n_release=len(r), max_abs_diff_eight=float(de.max()), max_abs_diff_all=float(d.max()), median_abs_diff_all=float(d.median())))
    rec['control_gate'] = dict(rows=rows, eight_max=e_max, all_max=a_max, eight_limit=GATE_EIGHT_MAX, all_limit=GATE_ALL_MAX)
    assert e_max <= GATE_EIGHT_MAX and a_max <= GATE_ALL_MAX, ('control (all KI nuclei) does not reproduce the S08 release: STOP and report', rows)
    return rel


def assemble(out, pins, gm, sel, rel, rec):
    eight = gm[gm.dataset_symbol.isin(MODEL_GENES)].set_index('stable_id')
    parts = []
    for sub in ('fit_control', 'fit_far'):
        idx = load_fit(out / sub); idx['fitdir'] = sub
        parts.append(idx)
    idx = pd.concat(parts, ignore_index=True)
    relg = rel.set_index(['stable_id', 'age_weeks'])
    seltot = sel.groupby(['rule', 'age']).agg(n_sel=('n_selected', 'sum'), n_tot=('n_nuclei', 'sum'))
    out_rows, full = [], []
    for row in idx.itertuples():
        tab = None
        if isinstance(row.table_file, str):
            tab = pd.read_csv(out / row.fitdir / 'tables' / row.table_file, sep='\t')
            f = tab.copy(); f.insert(0, 'variant', row.variant); f.insert(1, 'floor', row.floor); f.insert(2, 'fit', row.fit); f.insert(3, 'age', row.age)
            full.append(f)
        if row.age not in AGES:
            continue
        rule = {'control': 'all_KI_nuclei', 'wt90': 'wt_p90', 'top_third': 'top_third'}[row.variant]
        if row.variant == 'control':
            nk, share = int(sel[(sel.age == row.age) & (sel.rule == 'wt90')].n_nuclei.sum()), 1.0
        else:
            srule = 'wt90' if row.variant == 'wt90' else 'top_third'
            s = sel[(sel.age == row.age) & (sel.rule == srule)]
            if row.floor == FLOOR:
                s = s[s.n_selected >= FLOOR]
            nk = int(s.n_selected.sum()); share = nk / int(sel[(sel.age == row.age) & (sel.rule == srule)].n_nuclei.sum())
        for sid, g in eight.iterrows():
            rr = dict(gene=g.dataset_symbol, age=row.age, rule=rule, floor=row.floor, fit=row.fit, status=row.status, n_KI_mice=row.n_KI, n_WT_mice=row.n_WT, n_KI_nuclei=nk, share=share,
                      log2FC=np.nan, CI_low=np.nan, CI_high=np.nan)
            if tab is not None and sid in set(tab.stable_id):
                t = tab.set_index('stable_id').loc[sid]
                rr['log2FC'] = t.logFC
                if row.status == 'FITTED':
                    rr['CI_low'] = t['CI.L']; rr['CI_high'] = t['CI.R']
            rr['factor_mid'] = 2 ** rr['log2FC']; rr['factor_low'] = 2 ** rr['CI_low']; rr['factor_high'] = 2 ** rr['CI_high']
            if (sid, row.age) in relg.index:
                s8 = relg.loc[(sid, row.age)]
                rr.update(s08_perage_log2FC=s8.fallback_perage_log2FC, s08_perage_CI_low=s8.fallback_perage_CI_lower, s08_perage_CI_high=s8.fallback_perage_CI_upper,
                          s08_pooled_log2FC=s8.log2FC, s08_pooled_CI_low=s8.CI_lower, s08_pooled_CI_high=s8.CI_upper, s08_BH_q_decision=s8.BH_q)
            out_rows.append(rr)
    R = pd.DataFrame(out_rows)
    cols = ['gene', 'age', 'rule', 'log2FC', 'CI_low', 'CI_high', 'factor_mid', 'factor_low', 'factor_high', 'n_KI_nuclei', 'share', 'floor', 'fit', 'status', 'n_KI_mice', 'n_WT_mice',
            's08_perage_log2FC', 's08_perage_CI_low', 's08_perage_CI_high', 's08_pooled_log2FC', 's08_pooled_CI_low', 's08_pooled_CI_high', 's08_BH_q_decision']
    R = R[cols].sort_values(['gene', 'age', 'rule', 'fit', 'floor']).reset_index(drop=True)
    wtsv(R, out / 'FAR_SET_RATIOS.tsv', float_format='%.6g')
    wtsv(pd.concat(full, ignore_index=True), out / 'FULL_GENE_TABLES.tsv.gz', float_format='%.7g')
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', required=True, choices=['real', 'fixture'])
    ap.add_argument('--inputs'); ap.add_argument('--out', required=True)
    ap.add_argument('--rscript', required=True); ap.add_argument('--r-script', required=True)
    ap.add_argument('--r-sha256', required=True); ap.add_argument('--cleared-sha256')
    args = ap.parse_args()
    self_sha = sha(__file__)
    if args.mode == 'real':
        assert args.cleared_sha256 == self_sha, 'real mode needs --cleared-sha256 equal to this script'
    assert sha(args.r_script) == args.r_sha256, 'R script sha256 mismatch'
    out = Path(args.out); out.mkdir(parents=True, exist_ok=False)
    rec = dict(started_utc=datetime.now(timezone.utc).isoformat(), mode=args.mode, script_sha256=self_sha, r_script_sha256=args.r_sha256, disclaimer=DISCLAIMER,
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, platform=platform.platform(), floor=FLOOR,
               rules=dict(lane='unit vector of KI-WT difference of per-mouse mean scores over the separating lane PCs of the age; PC set re-derived from pc_check table and asserted',
                          wt90='KI nucleus score strictly greater than numpy-linear 90th percentile of the pooled WT nuclei of the same age',
                          top_third='ceil(n/3) highest-scoring KI nuclei pooled over the 3 KI mice of the age (stable order for ties)',
                          pseudobulk='sum of raw counts: far-set nuclei per KI mouse (18/24/30 W); all nuclei per WT mouse and for 5 W KI',
                          fits='24-column population matrix as in S08 (TMM + expression filter once, then 8-column pooled fit and per-age 3v3 fits with keep.lib.sizes); primary floor=1 (no exclusion), sensitivity floor=20'))
    pins = resolve_pins(args)
    rec['inputs'] = {k: dict(path=str(p), sha256=h) for k, (p, h) in pins.items()}
    cell, mice, sel, thr = stage_a(args, pins, out, rec)
    print('stage A done: SELECTION.tsv, OCCUPANCY.tsv, OCCUPANCY_SUMMARY.md written', flush=True)
    gm = stage_pseudobulk(pins, out, cell, mice, rec)
    run_r(args, pins, ['control'], out / 'pseudobulk', out / 'fit_control')
    rel = gate(out / 'fit_control', pins, gm, rec)
    print('control gate passed', rec['control_gate']['eight_max'], flush=True)
    run_r(args, pins, ['wt90', 'top_third'], out / 'pseudobulk', out / 'fit_far')
    R = assemble(out, pins, gm, sel, rel, rec)
    rec['completed_utc'] = datetime.now(timezone.utc).isoformat()
    wtext(out / 'RECORD.json', json.dumps(rec, indent=1, sort_keys=True, default=str) + '\n')
    lines = []
    for p in sorted(out.rglob('*')):
        if p.is_file() and p.name != 'SHA256SUMS':
            lines.append(f'{sha(p)}  {p.relative_to(out)}')
    wtext(out / 'SHA256SUMS', '\n'.join(lines) + '\n')
    print('done', out)


if __name__ == '__main__':
    main()

#%% debug
