"""P2 step (c): the five-gene fall in GSE75778 and GSE122099..

Primary values are the authors' Cuffdiff tables, cerebellum and inferior olive, mutant vs WT per line and age. Orientation is
set per table from its own FPKM columns (GSE75778: log2(fold_change) = log2(value_2 / value_1) with the mutant as sample_1,
so mutant vs WT = -log2(fold_change); GSE122099: Log2_Fold_Change = log2(mutant average / WT average)) and must hold for
>= 99.9 % of status-OK rows. Before any reading, the nine genes' mutant-vs-WT sign must equal the sign of the mean FPKM
difference (GSE75778: cuffnorm group means of the three replicates; GSE122099: the per-sample FPKMs in the same table): a
mismatch is a stop. The GSE75778 5-week labels A02, B05, D30 must correlate best with cuffnorm ATXN1[30Q], ATXN1[82Q] and
ATXN1[30Q]D776 (stop otherwise). Genes are matched by exact symbol (our gene_metadata symbols); in GSE122099 Kcnc3 exists only
as the merged locus "Kcnc3,Napsa", read both ways (merged locus used, flagged; Kcnc3 missing). Control genes: exact symbol,
one row, status OK, finite log2FC. CI: GSE75778 log2FC +- 1.96 |log2FC / test_stat|; GSE122099 has no test_stat, so |z| is
recovered from p as the normal quantile (approximate; its deviation from test_stat on GSE75778 is recorded).
Reading rule (pre-specified) applied mechanically to every comparison; roles as set at .
Read-only; writes a NEW real_v01/tests/ next to this script, exclusive create.

  python -I -B 02_replication/p2_bulk/p2_tests.py
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import binomtest, mannwhitneyu, norm

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
HERE = Path(__file__).resolve().parent
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/p2_tests'))
EXT = HERE / 'real_v01' / 'external_v02'
CONTROL = HERE / 'real_v01' / 'control'
PINS = {CONTROL / 'CONTROL_SET.tsv': '1968a1fdfbcb51371eb6d3cf1f725047d36fceeb790abd78fc412c5504bcc289',
        CONTROL / 'SHA256SUMS': '5f2043fee71304323223b47a9d76e4771ee3847a21ee6e8768fb79e19bcec021',
        EXT / 'SHA256SUMS': '127e70fadd1d457f38fdd7597ba5fc0531a37fa02faabbe83d426b82bdbc5305',
        ROOT / _release_paths.expand('${DATA_ROOT}/results/s08_pseudobulk_v01/gene_metadata.tsv'): '765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3'}
FIVE = ('Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g')
FOUR = ('Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2')
NINE = FIVE + FOUR
ORIENT_MIN, ORIENT_TOL, Z95, FALL_5W = 0.999, 1e-3, 1.959963984540054, -0.2
SIGN_GATE = 'threshold'    # 'threshold': a nine-gene sign mismatch stops only at |log2FC| >= SIGN_MIN_ABS; 'strict': any mismatch stops
SIGN_MIN_ABS, CONV_ABS, CONV_MIN = 0.1, 0.2, 0.999   # convention check: agreement over all OK genes with |log2FC| >= 0.2
F75 = 'GSE75778/files/GSE75778_'
F122 = 'GSE122099/files/GSE122099_'
CUFFNORM = F75 + 'cuffnorm_ATXN1.82Q_ATXN1.30Q.D776_WTFVB_genes.fpkm_tracking.txt.gz'
W5_LABELS = {'A02': 'ATXN1[30Q]', 'B05': 'ATXN1[82Q]', 'D30': 'ATXN1[30Q]D776'}
C75 = [  # file, sample_1 (mutant), sample_2 (WT), cuffnorm mutant group, age, line, role
    ('week5', 'B05', 'WT', 'ATXN1[82Q]', 5, 'ATXN1[82Q] Tg (B05)', '5 w: no fall yet'),
    ('week5', 'D30', 'WT', 'ATXN1[30Q]D776', 5, 'ATXN1[30Q]-D776 Tg', 'descriptive'),
    ('week5', 'A02', 'WT', 'ATXN1[30Q]', 5, 'ATXN1[30Q] Tg', 'descriptive'),
    ('week12', 'ATXN1[82Q]', 'FVB', 'ATXN1[82Q]', 12, 'ATXN1[82Q] Tg (B05)', 'secondary'),
    ('week12', 'ATXN1[30Q]D776', 'FVB', 'ATXN1[30Q]D776', 12, 'ATXN1[30Q]-D776 Tg', 'descriptive'),
    ('week28', 'ATXN1[82Q]', 'FVB', 'ATXN1[82Q]', 28, 'ATXN1[82Q] Tg (B05)', 'secondary'),
    ('week28', 'ATXN1[30Q]D776', 'FVB', 'ATXN1[30Q]D776', 28, 'ATXN1[30Q]-D776 Tg', 'descriptive')]
C122 = [  # file stem, tissue, age, line, role
    ('Atxn1_154Q_KI_Cerebellum_12weeks', 'cerebellum', 12, 'Atxn1 154Q/2Q KI', 'PRIMARY'),
    ('Atxn1_154Q_KI_Cerebellum_5weeks', 'cerebellum', 5, 'Atxn1 154Q/2Q KI', '5 w: no fall yet'),
    ('ATXN1_82Q_Tg_Cerebellum_12weeks', 'cerebellum', 12, 'ATXN1[82Q] Tg', 'secondary'),
    ('ATXN1_82Q_Tg_Cerebellum_5weeks', 'cerebellum', 5, 'ATXN1[82Q] Tg', '5 w: no fall yet'),
    ('Atxn1_154Q_KI_Inferior_Olive_12weeks', 'inferior olive', 12, 'Atxn1 154Q/2Q KI', 'descriptive'),
    ('Atxn1_154Q_KI_Inferior_Olive_5weeks', 'inferior olive', 5, 'Atxn1 154Q/2Q KI', 'descriptive'),
    ('ATXN1_82Q_Tg_Inferior_Olive_12weeks', 'inferior olive', 12, 'ATXN1[82Q] Tg', 'descriptive'),
    ('ATXN1_82Q_Tg_Inferior_Olive_5weeks', 'inferior olive', 5, 'ATXN1[82Q] Tg', 'descriptive')]
RULE_TEXT = ('P2 reading rule (pre-registered): At an age of 12 weeks or more, the fall “replicates” in a '
             'dataset if at least 4 of the 5 channel genes have log2FC < 0 and the one-sided sign test over the nine genes gives '
             'p < 0.05; it is “specific” if the median log2FC of the five channel genes is below the median of the control set '
             '(one-sided Wilcoxon rank-sum, p < 0.05). At 5 weeks the prediction from our data is no coordinated fall: fewer '
             'than 4 of 5 with log2FC < −0.2. Each dataset and age is reported against both parts of the rule. Caveats to state: '
             'bulk tissue, different constructs (82Q transgenic overexpression in Purkinje cells; 30Q-D776), different strains; '
             'direction and rank are the test, not effect size.')
NOTE = 'P2 tests: authors\' Cuffdiff values (GEO GSE75778, GSE122099), mutant vs WT; nine genes, control-set comparison and the plan\'s reading rule.'


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


def orientation(lfc, ratio):
    """+1 if lfc = ratio, -1 if lfc = -ratio, on >= 99.9 % of rows; otherwise a stop."""
    fwd, rev = float((np.abs(lfc - ratio) < ORIENT_TOL).mean()), float((np.abs(lfc + ratio) < ORIENT_TOL).mean())
    if fwd >= ORIENT_MIN:
        return 1, fwd, rev
    if rev >= ORIENT_MIN:
        return -1, fwd, rev
    raise ValueError(f'orientation not established (forward {fwd}, reverse {rev})')


def convention(frame, fpkm):
    """Share of status-OK genes with |log2FC| >= CONV_ABS whose sign equals the sign of the mean FPKM difference."""
    uniq = fpkm[~fpkm.symbol.duplicated(keep=False)].set_index('symbol')
    f = frame[(frame.status == 'OK') & np.isfinite(frame.lfc) & (frame.lfc.abs() >= CONV_ABS) & frame.symbol.isin(uniq.index)]
    if 'row' in fpkm:
        diff = (fpkm.mut - fpkm.wt).loc[f.index].to_numpy()
    else:
        diff = (uniq.mut - uniq.wt).reindex(f.symbol).to_numpy()
    agree = float((np.sign(f.lfc.to_numpy()) == np.sign(diff)).mean())
    if agree < CONV_MIN:
        raise ValueError(f'sign convention check failed: agreement {agree} over {len(f)} genes')
    return dict(convention_agreement=agree, convention_n=int(len(f)))


def load_75(ext):
    cn = pd.read_csv(ext / CUFFNORM, sep='\t')
    groups = {}
    for col in cn.columns:
        if col.endswith('_FPKM'):
            grp, rep = col[:-len('_FPKM')][:-1], col[:-len('_FPKM')][-1]
            groups.setdefault(grp, []).append(col)
    cn_mean = pd.DataFrame({g: cn[c].mean(axis=1) for g, c in groups.items()})
    cn_mean['symbol'] = cn.gene_short_name.astype(str)
    out, gates = [], {}
    tables = {w: pd.read_csv(ext / f'{F75}cuffdiff_{w}_ATXN1.82Q_ATXN1.30Q.D776_WTFVB_gene_exp.diff.txt.gz', sep='\t') for w in ('week5', 'week12', 'week28')}
    for w, mut, wt, grp, age, line, role in C75:
        d = tables[w]
        d = d[(d.sample_1 == mut) & (d.sample_2 == wt)].copy()
        L = pd.to_numeric(d.iloc[:, 9], errors='coerce')
        ok = (d.status == 'OK') & (d.value_1 > 0) & (d.value_2 > 0) & L.notna()
        sign, fwd, rev = orientation(L[ok].to_numpy(), np.log2(d.value_2[ok] / d.value_1[ok]).to_numpy())
        if sign != 1:
            raise ValueError(f'GSE75778 {w} {mut}: log2(fold_change) is not log2(value_2 / value_1)')
        wk = f'_wk{age:02d}'
        frame = pd.DataFrame(dict(symbol=d.gene.astype(str), status=d.status, lfc=-L, test_stat=-pd.to_numeric(d.test_stat, errors='coerce'),
                                  p=d.p_value, q=d.q_value, value_mut=d.value_1, value_wt=d.value_2)).reset_index(drop=True)
        fpkm = cn_mean[['symbol']].assign(mut=cn_mean[grp + wk], wt=cn_mean['FVB' + wk])
        out.append(dict(dataset='GSE75778', comparison=f'GSE75778_{w}_{mut}_vs_{wt}', tissue='cerebellum', age=age, line=line, role=role,
                        frame=frame, fpkm=fpkm, ci='test_stat'))
        gates[f'GSE75778_{w}_{mut}_vs_{wt}'] = dict(orientation_forward_fraction=fwd, orientation_reverse_fraction=rev, n_ok=int(ok.sum()),
                                                   **convention(frame, fpkm))
    # 5-week labels: each label's log2FC must correlate best with its assigned cuffnorm group
    d5 = tables['week5']
    corr = {}
    for lab, grp in W5_LABELS.items():
        x = d5[(d5.sample_1 == lab) & (d5.sample_2 == 'WT') & (d5.status == 'OK')]
        x = x.assign(lfc=-pd.to_numeric(x.iloc[:, 9], errors='coerce')).merge(cn_mean, left_on='gene', right_on='symbol')
        corr[lab] = {}
        for g in W5_LABELS.values():
            keep = (x[g + '_wk05'] > 0) & (x['FVB_wk05'] > 0) & np.isfinite(x.lfc)
            corr[lab][g] = float(np.corrcoef(x.lfc[keep], np.log2(x[g + '_wk05'][keep] / x['FVB_wk05'][keep]))[0, 1])
        if max(corr[lab], key=corr[lab].get) != grp:
            raise ValueError(f'5-week label {lab} does not correlate best with {grp}: {corr[lab]}')
    gates['GSE75778_week5_label_correlations'] = corr
    # recovered |z| from p against test_stat, on every GSE75778 comparison used
    ratios = []
    for c in out:
        f = c['frame']
        keep = (f.status == 'OK') & np.isfinite(f.test_stat) & (f.test_stat != 0) & (f.p > 0) & (f.p < 0.5)
        ratios.append(norm.isf(f.p[keep].to_numpy() / 2) / np.abs(f.test_stat[keep].to_numpy()))
    r = np.concatenate(ratios)
    gates['GSE75778_z_from_p_over_abs_test_stat'] = dict(n=int(r.size), q10=float(np.quantile(r, .1)), median=float(np.median(r)),
                                                         q90=float(np.quantile(r, .9)), p_floor=float(min(c['frame'].p[c['frame'].p > 0].min() for c in out)))
    return out, gates


def load_122(ext):
    out, gates = [], {}
    for stem, tissue, age, line, role in C122:
        d = pd.read_csv(ext / f'{F122}{stem}.txt.gz', sep='\t')
        cols = list(d.columns)
        gene, wt_reps, mut_reps, wt_avg, mut_avg = cols[1], cols[2:5], cols[5:8], cols[8], cols[9]
        if not (all(c.startswith('WT_') for c in wt_reps + [wt_avg]) and cols[10] == 'Test_Status' and cols[11] == 'Log2_Fold_Change' and cols[12] == 'p_value'):
            raise ValueError(f'unexpected GSE122099 columns: {stem} {cols}')
        num = {c: pd.to_numeric(d[c], errors='coerce') for c in wt_reps + mut_reps + [wt_avg, mut_avg, cols[11], cols[12], cols[13]]}
        L = num[cols[11]]
        ok = (d.Test_Status == 'OK') & (num[wt_avg] > 0) & (num[mut_avg] > 0) & L.notna()
        sign, fwd, rev = orientation(L[ok].to_numpy(), np.log2(num[mut_avg][ok] / num[wt_avg][ok]).to_numpy())
        frame = pd.DataFrame(dict(symbol=d[gene].astype(str), status=d.Test_Status, lfc=sign * L, test_stat=np.nan, p=num[cols[12]],
                                  q=num[cols[13]], value_mut=num[mut_avg], value_wt=num[wt_avg],
                                  lfc_raw=d[cols[11]].astype(str))).reset_index(drop=True)
        fpkm = pd.DataFrame(dict(symbol=frame.symbol, mut=pd.concat([num[c] for c in mut_reps], axis=1).mean(axis=1),
                                 wt=pd.concat([num[c] for c in wt_reps], axis=1).mean(axis=1), row=np.arange(len(frame))))
        out.append(dict(dataset='GSE122099', comparison=f'GSE122099_{stem}', tissue=tissue, age=age, line=line, role=role,
                        frame=frame, fpkm=fpkm, ci='from p (approximate)'))
        gates[f'GSE122099_{stem}'] = dict(orientation_sign=sign, orientation_forward_fraction=fwd, orientation_reverse_fraction=rev,
                                          n_ok=int(ok.sum()), non_numeric_log2fc=int(L.isna().sum()), q_column=cols[13], **convention(frame, fpkm))
    return out, gates


def find_rows(frame, symbol, allow_merged):
    hit = np.flatnonzero((frame.symbol == symbol).to_numpy())
    if len(hit) == 1:
        return int(hit[0]), ''
    if len(hit) > 1:
        return None, 'ambiguous: several rows'
    if allow_merged:
        merged = np.flatnonzero(frame.symbol.str.split(',').apply(lambda s: symbol in s).to_numpy())
        if len(merged) == 1:
            return int(merged[0]), f'merged locus {frame.symbol.iloc[merged[0]]}'
    return None, 'not found'


def nine_rows(c, allow_merged):
    rows = []
    f = c['frame']
    for g in NINE:
        i, flag = find_rows(f, g, allow_merged and g == 'Kcnc3')
        base = dict(dataset=c['dataset'], comparison=c['comparison'], tissue=c['tissue'], line=c['line'], age=c['age'], role=c['role'],
                    gene=g, gene_class='channel' if g in FIVE else 'other', row_symbol='', flag=flag)
        if i is None:
            rows.append(base)
            continue
        r = f.iloc[i]
        lfc = float(r.lfc)
        if c['ci'] == 'test_stat':
            z = abs(float(r.test_stat)) if np.isfinite(r.test_stat) and r.test_stat != 0 else np.nan
        else:
            z = float(norm.isf(r.p / 2)) if np.isfinite(r.p) and 0 < r.p < 1 else np.nan
        se = abs(lfc) / z if np.isfinite(lfc) and np.isfinite(z) and z > 0 else np.nan
        fp = c['fpkm']
        j = np.flatnonzero((fp.symbol == r.symbol).to_numpy()) if c['dataset'] == 'GSE75778' else np.array([i])
        fm, fw = (float(fp.mut.iloc[j[0]]), float(fp.wt.iloc[j[0]])) if len(j) == 1 else (np.nan, np.nan)
        base.update(row_symbol=r.symbol, status=r.status, log2fc_mut_vs_wt=lfc, test_stat=float(r.test_stat) if np.isfinite(r.test_stat) else np.nan,
                    abs_z_used=z, ci_lo=lfc - Z95 * se, ci_hi=lfc + Z95 * se, ci_method=c['ci'], p=float(r.p), q=float(r.q),
                    cuffdiff_value_mut=float(r.value_mut), cuffdiff_value_wt=float(r.value_wt), fpkm_mean_mut=fm, fpkm_mean_wt=fw,
                    n_fpkm_rows=int(len(j)))
        rows.append(base)
    return rows


def reading(c, nine, allow_merged):
    f = c['frame']
    found = nine[np.isfinite(nine.log2fc_mut_vs_wt.astype(float))] if 'log2fc_mut_vs_wt' in nine else nine.iloc[:0]
    five = found[found.gene_class == 'channel'].log2fc_mut_vs_wt.astype(float).to_numpy()
    allnine = found.log2fc_mut_vs_wt.astype(float).to_numpy()
    n_neg5, n_neg9 = int((five < 0).sum()), int((allnine < 0).sum())
    sign_p = float(binomtest(n_neg9, len(allnine), 0.5, alternative='greater').pvalue) if len(allnine) else np.nan
    counts = f.symbol.value_counts()
    ctrl = pd.read_csv(CONTROL / 'CONTROL_SET.tsv', sep='\t', comment='#')
    keep = f[f.symbol.isin(ctrl.symbol) & f.symbol.map(counts).eq(1) & (f.status == 'OK') & np.isfinite(f.lfc)]
    cl = keep.lfc.to_numpy(dtype=float)
    w_p = float(mannwhitneyu(five, cl, alternative='less').pvalue) if len(five) and len(cl) else np.nan
    flags = [f'{g}: {x}' for g, x in zip(nine.gene, nine.flag) if x]
    flags += [f'{g} status {s}' for g, s in zip(nine.gene, nine.status) if isinstance(s, str) and s != 'OK']
    flags += [f'{g} non-finite log2FC' for g, s, v in zip(nine.gene, nine.row_symbol, nine.log2fc_mut_vs_wt) if s and not np.isfinite(float(v))]
    r = dict(dataset=c['dataset'], comparison=c['comparison'], tissue=c['tissue'], line=c['line'], age=c['age'], role=c['role'],
             kcnc3=('merged locus Kcnc3,Napsa used' if allow_merged else 'missing') if c['dataset'] == 'GSE122099' else 'exact symbol',
             rule_part=('>= 12 w: replicates and specific' if c['age'] >= 12 else '5 w: no coordinated fall'),
             n_channel_found=len(five), n_channel_negative=n_neg5, n_nine_found=len(allnine), n_nine_negative=n_neg9, sign_test_p_one_sided=sign_p,
             replicates=bool(n_neg5 >= 4 and sign_p < 0.05), median_log2fc_channel=float(np.median(five)) if len(five) else np.nan,
             n_control=len(cl), median_log2fc_control=float(np.median(cl)) if len(cl) else np.nan, wilcoxon_p_one_sided=w_p,
             specific=bool(len(five) and len(cl) and np.median(five) < np.median(cl) and np.isfinite(w_p) and w_p < 0.05),
             n_channel_below_minus_0_2=int((five < FALL_5W).sum()),
             no_coordinated_fall=bool(int((five < FALL_5W).sum()) < 4), flags='; '.join(flags))
    return r, keep.assign(comparison=c['comparison'])[['comparison', 'symbol', 'status', 'lfc']]


def main():
    for path, want in PINS.items():
        if sha(path) != want:
            raise ValueError(f'input hash differs: {path}')
    listed = {name: digest for digest, name in (line.split('  ', 1) for line in (EXT / 'SHA256SUMS').read_text().splitlines())}
    used = [CUFFNORM] + [f'{F75}cuffdiff_{w}_ATXN1.82Q_ATXN1.30Q.D776_WTFVB_gene_exp.diff.txt.gz' for w in ('week5', 'week12', 'week28')] + [f'{F122}{c[0]}.txt.gz' for c in C122]
    for rel in used:
        if sha(EXT / rel) != listed[rel]:
            raise ValueError(f'external file differs from external_v02/SHA256SUMS: {rel}')
    comps75, g75 = load_75(EXT)
    comps122, g122 = load_122(EXT)
    comps = comps75 + comps122

    # sign assertion for the nine genes, before any reading
    nine_all, sign_check = [], []
    for c in comps:
        for allow in ((True, False) if c['dataset'] == 'GSE122099' else (False,)):
            rows = pd.DataFrame(nine_rows(c, allow))
            rows['kcnc3'] = ('merged' if allow else 'missing') if c['dataset'] == 'GSE122099' else 'exact'
            nine_all.append(rows)
    nine_tab = pd.concat(nine_all, ignore_index=True)
    for _, r in nine_tab[nine_tab.kcnc3 != 'missing'].iterrows():
        if not isinstance(r.get('row_symbol'), str) or not r.row_symbol:
            continue
        lfc, diff = float(r.log2fc_mut_vs_wt), float(r.fpkm_mean_mut - r.fpkm_mean_wt)
        testable = bool(np.isfinite(lfc))          # a non-finite log2FC is flagged and left out of the counts, not asserted
        ok = bool(np.isfinite(diff) and np.sign(lfc) == np.sign(diff)) if testable else None
        stops = bool(testable and not ok and (SIGN_GATE == 'strict' or abs(lfc) >= SIGN_MIN_ABS))
        sign_check.append(dict(comparison=r.comparison, gene=r.gene, row_symbol=r.row_symbol, log2fc_mut_vs_wt=lfc, fpkm_mean_mut=r.fpkm_mean_mut,
                               fpkm_mean_wt=r.fpkm_mean_wt, n_fpkm_rows=r.n_fpkm_rows, testable=testable, sign_matches=ok, stops=stops,
                               flag='' if ok or not testable or stops else f'sign mismatch at |log2FC| < {SIGN_MIN_ABS} (flagged, not a stop)'))
    sign_check = pd.DataFrame(sign_check)
    bad = sign_check[sign_check.stops]
    if len(bad):
        print(bad.to_string(index=False))
        raise ValueError('FPKM sign assertion failed for the rows above: STOP before any reading')

    readings, control_rows = [], []
    for c in comps:
        for allow in ((True, False) if c['dataset'] == 'GSE122099' else (False,)):
            key = ('merged' if allow else 'missing') if c['dataset'] == 'GSE122099' else 'exact'
            nine = nine_tab[(nine_tab.comparison == c['comparison']) & (nine_tab.kcnc3 == key)].reset_index(drop=True)
            r, ctrl = reading(c, nine, allow)
            readings.append(r)
            if key != 'missing':
                control_rows.append(ctrl)
    read = pd.DataFrame(readings)
    ctrl_tab = pd.concat(control_rows, ignore_index=True)

    lines = ['# P2 five-gene fall in GSE75778 and GSE122099', '', NOTE, '',
             '| comparison | role | Kcnc3 | channel neg / found | nine neg / found | sign p | replicates | median channel | median control (n) | Wilcoxon p | specific | 5 w: n < -0.2 | no fall (5 w) |',
             '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    yn = lambda v: 'yes' if v else 'no'
    for r in read.itertuples():
        star = '' if r.age >= 12 else ' (n/a at 5 w)'
        lines.append(f'| {r.comparison} | {r.role} | {r.kcnc3} | {r.n_channel_negative}/{r.n_channel_found} | {r.n_nine_negative}/{r.n_nine_found} | '
                     f'{r.sign_test_p_one_sided:.4g} | {yn(r.replicates)}{star} | {r.median_log2fc_channel:.3f} | {r.median_log2fc_control:.3f} ({r.n_control}) | '
                     f'{r.wilcoxon_p_one_sided:.4g} | {yn(r.specific)}{star} | {r.n_channel_below_minus_0_2} | {yn(r.no_coordinated_fall)}{"" if r.age < 12 else " (n/a at >= 12 w)"} |')
    zr = g75['GSE75778_z_from_p_over_abs_test_stat']
    lines += ['', f'CI note: GSE122099 has no test_stat; |z| is recovered from p (normal quantile) and the CI is approximate. On GSE75778 the '
                  f'recovered |z| is {zr["median"]:.2f} x |test_stat| (median; q10-q90 {zr["q10"]:.2f}-{zr["q90"]:.2f}; Cuffdiff p floor '
                  f'{zr["p_floor"]:.0e}), so these CIs are about that much narrower than the GSE75778 ones and are not comparable.', '',
              'Caveats: bulk tissue; different constructs (82Q transgenic overexpression in Purkinje cells; 30Q-D776); different strains '
              '(FVB for the transgenics, C57BL/6J for the knock-in); direction and rank are the test, not effect size.', '', RULE_TEXT, '']

    OUT.mkdir(parents=False, exist_ok=False)
    wtsv(nine_tab, 'NINE_GENES.tsv')
    wtsv(read, 'READING.tsv')
    wtsv(sign_check, 'SIGN_ASSERTION.tsv')
    wtsv(ctrl_tab, 'CONTROL_LOG2FC.tsv')
    with open(OUT / 'SUMMARY.md', 'x') as f:
        f.write('\n'.join(lines))
    rec = dict(note=NOTE, task='P2 step (c)', script_sha256=sha(Path(__file__)), inputs={str(p.relative_to(ROOT)): h for p, h in PINS.items()},
               external_files_used={rel: listed[rel] for rel in used}, gates={**g75, **g122}, rule_text=RULE_TEXT,
               sign_gate=dict(mode=SIGN_GATE, min_abs_log2fc=SIGN_MIN_ABS, convention_abs=CONV_ABS, convention_min=CONV_MIN,
                              flagged_small_mismatches=sign_check[sign_check.flag != ''][['comparison', 'gene', 'log2fc_mut_vs_wt']].to_dict('records')),
               choices=dict(kcnc3_gse122099='read both ways: merged locus Kcnc3,Napsa used (flagged) and Kcnc3 missing',
                            ci='GSE75778 +-1.96 |log2FC/test_stat|; GSE122099 |z| from p, approximate', control='exact symbol, one row, status OK, finite log2FC',
                            sign_test='binomtest(n negative, n found, 0.5, greater)', wilcoxon="mannwhitneyu(channel, control, alternative='less')",
                            roles='PRIMARY 154Q/2Q cerebellum 12 w; secondary 82Q cerebellum 12 w (both series) and 28 w; 5 w no-fall test; IO, 30Q-D776, 30Q descriptive'),
               python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
               finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)) + '\n')
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.name}\n')
    print('\n'.join(lines[:4 + len(read) + 2]))
    print(json.dumps({p.name: sha(p) for p in files}, indent=1))


if __name__ == '__main__':
    main()

#%% debug
