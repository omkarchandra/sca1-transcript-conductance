"""Phase 4 Track 2 , backlog B3 (pre-specified): convert the S08 limma release into per-gene KI/WT expression
ratios for the tested M1-M5 Purkinje genes per age, then into conductance scale factors under the spec's mapping rule.

RULE (stated, with its source): TRACK2_SPEC_v01 section 3.3 "Scaling rule" (unchanged by spec v02) defines the exact primary
conductance factor from per-mouse CPM: sum of the tested pore-gene CPM in a mouse divided by the corresponding WT reference
sum. The S08 release has per-gene `log2FC_decision` (limma/voom normalized logCPM contrast) and `WT_raw_sum` / `KI_raw_sum`
(unnormalized counts), but not the per-sample effective library sizes or WT pseudobulk CPM needed to weight the family sum.
Therefore this file outputs (i) per-gene normalized KI/WT effect ratios = 2^log2FC_decision with exponentiated decision-CI
bounds; (ii) a separately labelled **raw-count family-ratio proxy** = sum KI_raw_sum / sum WT_raw_sum over tested pore genes;
and (iii) `CPM_PRIMARY_FACTOR_UNAVAILABLE` rather than claiming that proxy is the spec's exact factor. This proxy is
DESCRIPTIVE ONLY and must not be used as a worker scale setting; the exact primary factor is computed from the pinned Step 8
per-mouse pseudobulk CPM when that acceptance run is . This distinction is essential: weighting limma ratios by raw
counts would mix normalized estimates with unnormalized weights.
Genes with a non-empty `failure` or without a numeric log2FC_decision are EXCLUDED from the factor and listed (NO_ESTIMATE);
genes missing an age are listed (MISSING_AGE); ratios are NOT clipped and every ratio outside [0.1, 10] is counted and listed
(section 3.3 verbatim). G_c = pore-forming genes of the conductance table (p4_track2_gene_tables.FAMILIES) intersected with
the module's TESTED symbols in FREEZE.json (symbol_table minus absent_or_filtered of the Purkinje entry; definitions only).
ADAPTATION, stated plainly: section 3.3's PRIMARY rule is per-mouse pseudobulk (Step 8 CPM) with a leave-one-out WT reference;
the release is aggregate (3 KI vs 3 WT per age), so this converter produces the RELEASE-LEVEL gene-effect view and a raw-count
family proxy for documentation, but cannot produce the exact CPM primary factor from the release alone. Per-mouse factors
keep using the Step 8 pseudobulk when those runs are . Auxiliary subunits
(Atp1b*, Fxyd*) are listed and NOT scaled (section 3.3 sensitivity S3). b_Ca is sensitivity-only; g_ERG is inert (gbar 0).
NO REAL RUN before review: this file was delivered with fixtures and a baselined mutation pass only; opening the pinned
inputs at run time re-verifies their sha256 first and refuses anything else.

Usage: python3 -I ${REPO_ROOT}/05_model_forrest/p4_t2_scale_factors_v01.py <release .tsv.gz> <FREEZE .json> <new output dir>
Writes GENE_RATIOS.tsv, CONDUCTANCE_FACTORS.tsv, AUXILIARY.tsv and RESULT.json (all exclusive-create) into the new output dir.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys

SCRIPTS=Path(__file__).resolve().parent
import types

RELEASE_SHA='40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298'
FREEZE_SHA='aaae299efb5bceb4ba0363ee6962c77f1b3cdf11479aec0be1e600efd8f8a2a5'
TABLES_PATH=SCRIPTS/'p4_track2_gene_tables.py'
TABLES_SHA='837dfe286bcb945252f3b16e83198e1917b3940b0ed2ecab0645b844ea86cc3d'
PINS={'tables':TABLES_SHA,'release':RELEASE_SHA,'freeze':FREEZE_SHA}
AGES=(5,18,24,30)
POPULATION='purkinje'
tables=None

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_tables(pins):
    """Verify the imported family-table bytes before executing them or opening release/FREEZE inputs."""
    data=TABLES_PATH.read_bytes()
    got=hashlib.sha256(data).hexdigest()
    if got!=pins['tables']:raise ValueError(f'gene table module differs from its pin: {got}')
    name='p4_track2_gene_tables'
    module=sys.modules.get(name)
    if module is None:
        module=types.ModuleType(name);module.__file__=str(TABLES_PATH)
        exec(compile(data,str(TABLES_PATH),'exec'),module.__dict__)
        sys.modules[name]=module
    return module,got

def get_tables():
    global tables
    if tables is None:tables,_=load_tables(PINS)
    return tables

def num(text):
    try:
        value=float(text)
        return value if math.isfinite(value) else None
    except (TypeError,ValueError):return None

def read_release(path):
    rows=[];seen=set()
    with gzip.open(path,'rt') as stream:
        header=stream.readline().rstrip('\n').split('\t')
        required=('dataset_symbol','population','age_weeks','log2FC_decision','CI_lower_decision','CI_upper_decision','WT_raw_sum','KI_raw_sum','n_WT','failure','decision_source','stable_id')
        for name in required:
            if name not in header:raise ValueError(f'release column missing: {name}')
        for line in stream:
            parts=line.rstrip('\n').split('\t')
            if len(parts)!=len(header):raise ValueError('ragged release row')
            row=dict(zip(header,parts))
            key=(row['dataset_symbol'],row['population'],row['age_weeks'])
            if key in seen:raise ValueError(f'duplicate symbol/population/age release row: {key}')
            seen.add(key);rows.append(row)
    return rows

def tested_symbols(freeze):
    """Per module: the Purkinje TESTED symbol set (definitions only): symbol_table minus absent_or_filtered."""
    out={}
    for module,block in freeze.items():
        purkinje=block.get('per_population',{}).get('Purkinje cell')
        if purkinje is None:continue
        resolved={entry['symbol'] for entry in block.get('symbol_table',[])}
        excluded={entry['symbol'] if isinstance(entry,dict) else entry for entry in (purkinje.get('absent_or_filtered') or [])}
        out[module]=resolved-excluded
    return out

def csv_members(value):
    """Parse a sorted comma field, treating the explicit '-' placeholder as empty."""
    return set(value.split(',')) - {'','-'}

def convert(rows,freeze):
    table=get_tables()
    tested=tested_symbols(freeze)
    for conductance,_,module,_,_,_ in table.FAMILIES:
        if module not in tested:raise ValueError(f'module of the conductance table missing in the freeze: {module}')
    wanted={}
    for conductance,_,module,pore,aux,_ in table.FAMILIES:
        wanted[conductance]=(module,set(pore),set(aux))
    release_by_gene={}
    for row in rows:
        if row['population']!=POPULATION:continue
        age=num(row['age_weeks'])
        if age is None or not age.is_integer() or int(age) not in AGES:raise ValueError(f'age outside the frozen four: {row["age_weeks"]}')
        release_by_gene.setdefault(row['dataset_symbol'],{})[int(age)]=row
    gene_rows=[];problems=[]
    for conductance,(module,pore,aux) in wanted.items():
        for gene in sorted(pore|aux):
            per_age=release_by_gene.get(gene,{})
            for age in AGES:
                row=per_age.get(age)
                if row is None:
                    if gene in pore:problems.append({'conductance':conductance,'gene':gene,'age':age,'problem':'MISSING_AGE'})
                    continue
                estimate=num(row['log2FC_decision']);low=num(row['CI_lower_decision']);high=num(row['CI_upper_decision'])
                failed=bool(row['failure'].strip())
                wt_raw=num(row['WT_raw_sum']);ki_raw=num(row['KI_raw_sum'])
                if wt_raw is not None and wt_raw<0:raise ValueError(f'negative WT raw sum for {gene} at {age}w')
                if ki_raw is not None and ki_raw<0:raise ValueError(f'negative KI raw sum for {gene} at {age}w')
                entry={'conductance':conductance,'module':module,'gene':gene,'role':'pore' if gene in pore else 'auxiliary',
                       'age_weeks':age,'tested':gene in tested[module],'decision_source':row['decision_source'],
                       'n_WT':row['n_WT'],'WT_raw_sum':wt_raw,'KI_raw_sum':ki_raw,
                       'log2FC_decision':estimate,'ratio':None,'ratio_ci_low':None,'ratio_ci_high':None,
                       'raw_count_ratio':ki_raw/wt_raw if wt_raw is not None and wt_raw>0 and ki_raw is not None else None,
                       'problem':None}
                if failed or estimate is None:
                    entry['problem']='NO_ESTIMATE' if gene in tested[module] and gene in pore else 'EXCLUDED_UNTESTED_OR_AUX'
                else:
                    try:
                        entry['ratio']=2.0**estimate
                        entry['ratio_ci_low']=2.0**low if low is not None else None
                        entry['ratio_ci_high']=2.0**high if high is not None else None
                    except OverflowError as error:raise ValueError(f'ratio overflow for {gene} at {age}w') from error
                gene_rows.append(entry)
                if entry['problem'] and 'EXCLUDED' not in (entry['problem'] or ''):problems.append({'conductance':conductance,'gene':gene,'age':age,'problem':entry['problem']})
    factor_rows=[]
    for conductance,(module,pore,aux) in wanted.items():
        used=[r for r in gene_rows if r['conductance']==conductance and r['role']=='pore' and r['tested'] and r['problem'] is None]
        genes_tested=sorted(gene for gene in pore if gene in tested[module])
        factor={'conductance':conductance,'module':module,'genes_tested':','.join(genes_tested) or '-','n_genes_tested':len(genes_tested)}
        for age in AGES:
            at_age=[r for r in used if r['age_weeks']==age]
            excluded=sorted({p['gene'] for p in problems if p['conductance']==conductance and p['age']==age})
            column=f'raw_count_proxy_{age}w'
            factor[f'missing_raw_counts_{age}w']='-'
            if not genes_tested:
                factor[column]=None;factor[f'proxy_status_{age}w']='NO_TESTED_GENE'
                factor[f'cpm_primary_factor_{age}w']=None;factor[f'cpm_status_{age}w']='NO_TESTED_GENE'
                factor[f'included_genes_{age}w']='-'
            elif not at_age:
                factor[column]=None;factor[f'proxy_status_{age}w']='NO_ESTIMATE_AT_AGE'
                factor[f'cpm_primary_factor_{age}w']=None;factor[f'cpm_status_{age}w']='NO_ESTIMATE_AT_AGE'
                factor[f'included_genes_{age}w']='-'
            else:
                missing_raw=sorted(r['gene'] for r in at_age if r['WT_raw_sum'] is None or r['KI_raw_sum'] is None)
                factor[f'missing_raw_counts_{age}w']=','.join(missing_raw) or '-'
                if missing_raw:
                    factor[column]=None;factor[f'proxy_status_{age}w']='MISSING_RAW_COUNTS'
                    factor[f'included_genes_{age}w']='-'
                else:
                    total_wt=sum(r['WT_raw_sum'] for r in at_age)
                    total_ki=sum(r['KI_raw_sum'] for r in at_age)
                    factor[f'included_genes_{age}w']=','.join(sorted(r['gene'] for r in at_age)) or '-'
                    if total_wt<=0:factor[column]=None;factor[f'proxy_status_{age}w']='ZERO_REFERENCE'
                    else:
                        factor[column]=total_ki/total_wt;factor[f'proxy_status_{age}w']='OK_RAW_COUNTS_ONLY'
                factor[f'cpm_primary_factor_{age}w']=None;factor[f'cpm_status_{age}w']='CPM_PRIMARY_FACTOR_UNAVAILABLE_FROM_RELEASE'
            factor[f'excluded_{age}w']=','.join(excluded) or '-'
        factor['ratios_outside_0.1_10']=[r['gene']+f'@{r["age_weeks"]}w' for r in used if r['ratio'] is not None and not 0.1<=r['ratio']<=10]
        factor['note']='sensitivity only' if conductance=='b_Ca' else ('inert (gbar 0)' if conductance=='g_ERG' else '')
        factor_rows.append(factor)
    return gene_rows,factor_rows,problems

def write_tsv(path,rows,columns):
    with Path(path).open('x') as stream:
        stream.write('\t'.join(columns)+'\n')
        for row in rows:stream.write('\t'.join('' if row.get(c) is None else (str(row[c]) if not isinstance(row.get(c),list) else ','.join(row[c])) for c in columns)+'\n')

def main(argv=None):
    args=list(sys.argv[1:] if argv is None else argv)
    if len(args)!=3:raise SystemExit('expected exactly: <release .tsv.gz> <FREEZE .json> <new output dir>')
    release,freeze_path,out=[Path(a) for a in args]
    global tables
    tables,tables_sha=load_tables(PINS)
    out.mkdir(parents=True,exist_ok=False)
    pins=dict(PINS)
    checks=[{'check':'gene table module sha256 matches the pin','passed':True,'got':tables_sha}]
    for name,path in (('release',release),('freeze',freeze_path)):
        got=sha(path);checks.append({'check':f'{name} sha256 matches the pin','passed':got==pins[name],'got':got})
        if got!=pins[name]:raise ValueError(f'{name} differs from its pin: {got}')
    rows=read_release(release)
    freeze=json.loads(freeze_path.read_text())
    gene_rows,factor_rows,problems=convert(rows,freeze)
    recomputed=0
    for factor in factor_rows:
        for age in AGES:
            value=factor[f'raw_count_proxy_{age}w']
            if value is None:continue
            at_age=[r for r in gene_rows if r['conductance']==factor['conductance'] and r['role']=='pore' and r['tested'] and r['problem'] is None and r['age_weeks']==age]
            total_wt=sum(r['WT_raw_sum'] or 0.0 for r in at_age);total_ki=sum(r['KI_raw_sum'] or 0.0 for r in at_age)
            if total_wt<=0 or abs(total_ki/total_wt-value)>1e-12*max(1.0,abs(value)):raise ValueError('self-check: raw-count proxy recomputation differs')
            recomputed+=1
    checks.append({'check':'every non-null raw-count proxy recomputes from the per-gene WT/KI raw sums','passed':True,'n':recomputed})
    checks.append({'check':'every included family gene belongs to that conductance\'s tested pore set','passed':all(csv_members(factor.get(f'included_genes_{age}w','-'))<=csv_members(factor['genes_tested']) for factor in factor_rows for age in AGES)})
    checks.append({'check':'every exact CPM primary factor is explicitly unavailable from this release','passed':all(factor.get(f'cpm_primary_factor_{age}w') is None for factor in factor_rows for age in AGES)})
    checks.append({'check':'every CPM status matches a defined unavailable or empty-gene state','passed':all(factor.get(f'cpm_status_{age}w') in ('CPM_PRIMARY_FACTOR_UNAVAILABLE_FROM_RELEASE','NO_TESTED_GENE','NO_ESTIMATE_AT_AGE') for factor in factor_rows for age in AGES)})
    gene_columns=('conductance','module','gene','role','tested','age_weeks','decision_source','n_WT','WT_raw_sum','KI_raw_sum','log2FC_decision','ratio','ratio_ci_low','ratio_ci_high','raw_count_ratio','problem')
    factor_columns=('conductance','module','genes_tested','n_genes_tested')+tuple(v for age in AGES for v in (f'raw_count_proxy_{age}w',f'proxy_status_{age}w',f'cpm_primary_factor_{age}w',f'cpm_status_{age}w',f'included_genes_{age}w',f'missing_raw_counts_{age}w',f'excluded_{age}w'))+('ratios_outside_0.1_10','note')
    write_tsv(out/'GENE_RATIOS.tsv',gene_rows,gene_columns)
    write_tsv(out/'CONDUCTANCE_FACTORS.tsv',factor_rows,tuple(factor_columns))
    aux=[r for r in gene_rows if r['role']=='auxiliary']
    write_tsv(out/'AUXILIARY.tsv',aux,gene_columns)
    result={'script_sha256':sha(__file__),'tables_sha256':tables_sha,'release_sha256':sha(release),'freeze_sha256':sha(freeze_path),
            'rule':'TRACK2_SPEC_v01 section 3.3 primary is the WT-referenced sum of per-mouse CPM across tested pore genes; this release lacks per-mouse CPM and normalized WT weights, so only raw-count family proxies are emitted and exact CPM scale factors are explicitly unavailable; per-gene effect ratios are 2^log2FC_decision',
            'ages':AGES,'n_gene_rows':len(gene_rows),'n_factors':len(factor_rows),'n_problems':len(problems),'problems':problems,
            'n_checks':len(checks),'n_pass':sum(c['passed'] for c in checks),'all_pass':all(c['passed'] for c in checks),'checks':checks}
    with (out/'RESULT.json').open('x') as stream:json.dump(result,stream,indent=1)
    print(result['n_pass'],'/',result['n_checks'],'all_pass',result['all_pass'],'| factors',len(factor_rows),'| problems',len(problems))
    sys.exit(0 if result['all_pass'] else 1)

if __name__=='__main__':main()

#%% debug
