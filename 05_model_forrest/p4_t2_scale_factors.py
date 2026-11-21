"""B3 v02 release-level factors .

This is a new, append-only derivation. It leaves the B3 v01 files untouched.
The normalised effect-weighted family factor for conductance c and age a is

    F_norm(c,a) = sum_g WT_raw_sum(g,a) * 2**log2FC_decision(g,a)
                  / sum_g WT_raw_sum(g,a)

over tested pore genes with usable estimates. WT_raw_sum is an unnormalised
count weight, so F_norm is explicitly a release-level normalized-effect /
raw-weight proxy sensitivity, not a normalized-CPM family factor and not the
spec's exact per-mouse CPM primary factor. The unweighted geometric mean is
reported separately as the predeclared S1 sensitivity. No family confidence
interval is claimed. Neither aggregate is a worker scale setting.

The v01 raw-count proxy columns are retained byte-for-byte. Exact CPM primary
factors remain unavailable. b_Ca is sensitivity-only: tested Calb1/Pvalb
rows absent from the S08 release are sourced from the pinned complete S08
table and labelled with their file and family. Nothing here loads or runs the
Purkinje model.

Usage (real data only after independent byte review and a separate run
approval):
  python -I ${REPO_ROOT}/05_model_forrest/p4_t2_scale_factors.py     <S08 release.tsv.gz> <FREEZE.json> <complete S08 results.tsv.gz>     <v01 CONDUCTANCE_FACTORS.tsv> <new output dir>

The planned real output directory is
${DATA_ROOT}/results/p4_t2_scale_factors_v02. All output files are exclusive-
create; do not rerun into an existing directory.
"""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import sys

SCRIPTS = Path(__file__).resolve().parent
ROOT = _ReleasePath(_release_paths.expand('${REPO_ROOT}'))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
import p4_t2_scale_factors_v01 as v1

AGES = (5, 18, 24, 30)
V1_SCRIPT_SHA = '6af5783be1e37bc85223de76ed0a85ab1c24c687b0b94e745c86936e21e4daf7'
TABLES_SHA = '837dfe286bcb945252f3b16e83198e1917b3940b0ed2ecab0645b844ea86cc3d'
RELEASE_SHA = '40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298'
FREEZE_SHA = 'aaae299efb5bceb4ba0363ee6962c77f1b3cdf11479aec0be1e600efd8f8a2a5'
COMPLETE_SHA = 'e66e4ae4fdd4a2c101df45c893ca25308e9c650e66255990095effdf4655c44a'
V1_FACTORS_SHA = '09ac0183c64f0a6df2e2c125c73f7aeb5295d6af6e384d8f432c73a87e55829c'

RELEASE_REL = _release_paths.expand('${REPO_ROOT}/data/snrnaseq/S08_S11_release_purkinje_bergmann.tsv.gz')
FREEZE_REL = _release_paths.expand('${DATA_ROOT}/s08_modules_freeze_v03/FREEZE.json')
COMPLETE_REL = _release_paths.expand('${DATA_ROOT}/s08_limma_real_results_v02_derived/expression_results_complete.tsv.gz')
V1_FACTORS_REL = _release_paths.expand('${DATA_ROOT}/results/p4_t2_scale_factors_v01/CONDUCTANCE_FACTORS.tsv')
PLANNED_REAL_OUT = _release_paths.expand('${DATA_ROOT}/results/p4_t2_scale_factors_v02')

PINNED = {
    'v1_converter': V1_SCRIPT_SHA,
    'tables': TABLES_SHA,
    'release': RELEASE_SHA,
    'freeze': FREEZE_SHA,
    'complete_results': COMPLETE_SHA,
    'v1_factors': V1_FACTORS_SHA,
}
COMPLETE_REQUIRED = (
    'stable_id', 'dataset_symbol', 'population', 'age_weeks', 'family',
    'log2FC_decision', 'CI_lower_decision', 'CI_upper_decision',
    'WT_raw_sum', 'KI_raw_sum', 'n_WT', 'failure', 'decision_source',
)
V1_GENE_COLUMNS = (
    'conductance', 'module', 'gene', 'role', 'tested', 'age_weeks',
    'decision_source', 'n_WT', 'WT_raw_sum', 'KI_raw_sum',
    'log2FC_decision', 'ratio', 'ratio_ci_low', 'ratio_ci_high',
    'raw_count_ratio', 'problem',
)
GENE_COLUMNS = V1_GENE_COLUMNS + ('stable_id', 'source_file', 'source_family')
NORM_SUFFIXES = (
    'f_norm', 'f_norm_status', 'f_norm_weight_sum',
    'f_norm_included_genes', 'f_norm_excluded_genes', 'f_norm_source_families',
    'geometric_mean_sensitivity', 'geometric_mean_status',
    'geometric_mean_included_genes', 'geometric_mean_excluded_genes',
)
BUFFER_COLUMNS = (
    'conductance', 'module', 'gene', 'stable_id', 'age_weeks',
    'log2FC_decision', 'ratio', 'ratio_ci_low', 'ratio_ci_high',
    'WT_raw_sum', 'KI_raw_sum', 'decision_source', 'tested', 'problem',
    'source_file', 'source_family',
)


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def age_int(value):
    parsed = number(value)
    if parsed is None or not parsed.is_integer() or int(parsed) not in AGES:
        raise ValueError(f'age outside frozen four: {value}')
    return int(parsed)


def row_text(value):
    if value is None:
        return ''
    if isinstance(value, (list, tuple, set)):
        return ','.join(str(v) for v in value)
    if isinstance(value, bool):
        return 'True' if value else 'False'
    return str(value)


def write_tsv_exclusive(path, rows, columns):
    with Path(path).open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='raise')
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row_text(row.get(key)) for key in columns})


def read_tsv(path):
    with Path(path).open('r', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames:
            raise ValueError(f'missing TSV header: {path}')
        rows = list(reader)
    return list(reader.fieldnames), rows


def registry_population(release_row):
    name = release_row.get('population_registry_name', '').strip()
    if name:
        return name
    return {'purkinje': 'Purkinje cell', 'bergmann': 'Bergmann glia'}.get(
        release_row.get('population', ''), release_row.get('population', ''))


def release_row_keys(release_rows):
    keys = set()
    for row in release_rows:
        if row.get('family') != 'core_primary':
            raise ValueError('S08 release contains a non-core_primary row')
        key = (row['stable_id'], registry_population(row), age_int(row['age_weeks']))
        if key in keys:
            raise ValueError(f'duplicate release identity: {key}')
        keys.add(key)
    return keys


def index_complete_rows(complete_rows, release_rows, supplemental_symbols):
    """Index only release matches and b_Ca secondary rows; all records are caller-supplied."""
    wanted_core = release_row_keys(release_rows)
    wanted_core_prefixes = {(stable_id, population) for stable_id, population, _ in wanted_core}
    wanted_secondary = {(symbol, 'Purkinje cell', age)
                        for symbol in supplemental_symbols for age in AGES}
    core = {}
    secondary = {}
    for source in complete_rows:
        family = source.get('family')
        if family not in ('core_primary', 'whole_transcriptome_secondary'):
            continue
        population = source.get('population')
        if family == 'core_primary' and (source.get('stable_id'), population) not in wanted_core_prefixes:
            continue
        if family == 'whole_transcriptome_secondary' and (
                population != 'Purkinje cell' or source.get('dataset_symbol') not in supplemental_symbols):
            continue
        age = age_int(source.get('age_weeks'))
        if family == 'core_primary':
            key = (source['stable_id'], population, age)
            if key in wanted_core:
                if key in core:
                    raise ValueError(f'duplicate complete core identity: {key}')
                core[key] = source
        elif population == 'Purkinje cell' and source.get('dataset_symbol') in supplemental_symbols:
            key = (source['dataset_symbol'], population, age)
            if key in wanted_secondary:
                if key in secondary:
                    raise ValueError(f'duplicate complete secondary identity: {key}')
                secondary[key] = source
    missing = sorted(wanted_core - set(core))
    if missing:
        raise ValueError(f'complete S08 table lacks release core rows: {missing[:5]}')
    return core, secondary


def read_complete_index(path, release_rows, supplemental_symbols):
    """Stream the pinned complete table; retain only rows needed for B3 checks/outputs."""
    wanted_core = release_row_keys(release_rows)
    wanted_core_prefixes = {(stable_id, population) for stable_id, population, _ in wanted_core}
    wanted_secondary = {(symbol, 'Purkinje cell', age)
                        for symbol in supplemental_symbols for age in AGES}
    core = {}
    secondary = {}
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('complete S08 table has no header')
        missing = sorted(set(COMPLETE_REQUIRED) - set(reader.fieldnames))
        if missing:
            raise ValueError(f'complete S08 header missing: {missing}')
        for source in reader:
            family = source['family']
            if family not in ('core_primary', 'whole_transcriptome_secondary'):
                continue
            population = source['population']
            if family == 'core_primary':
                if (source['stable_id'], population) not in wanted_core_prefixes:
                    continue
                age = age_int(source['age_weeks'])
                key = (source['stable_id'], population, age)
                if key in wanted_core:
                    if key in core:
                        raise ValueError(f'duplicate complete core identity: {key}')
                    core[key] = source
            elif population == 'Purkinje cell' and source['dataset_symbol'] in supplemental_symbols:
                age = age_int(source['age_weeks'])
                key = (source['dataset_symbol'], population, age)
                if key in wanted_secondary:
                    if key in secondary:
                        raise ValueError(f'duplicate complete secondary identity: {key}')
                    secondary[key] = source
    missing_core = sorted(wanted_core - set(core))
    if missing_core:
        raise ValueError(f'complete S08 table lacks release core rows: {missing_core[:5]}')
    return core, secondary


def check_release_complete_identity(release_rows, complete_core):
    """Require exact agreement on the effects and sums used from every release row."""
    compared = 0
    fields = ('stable_id', 'dataset_symbol', 'n_WT', 'failure', 'decision_source',
              'log2FC_decision', 'CI_lower_decision', 'CI_upper_decision',
              'WT_raw_sum', 'KI_raw_sum')
    for release in release_rows:
        key = (release['stable_id'], registry_population(release), age_int(release['age_weeks']))
        source = complete_core.get(key)
        if source is None:
            raise ValueError(f'missing complete-table counterpart for release row: {key}')
        for field in fields:
            if release.get(field, '') != source.get(field, ''):
                raise ValueError(f'release/complete mismatch for {key} field {field}: '
                                 f'{release.get(field)!r} != {source.get(field)!r}')
        if source.get('family') != 'core_primary':
            raise ValueError(f'release counterpart is not core_primary: {key}')
        compared += 1
    return compared


def read_v1_factor_file(path):
    columns, rows = read_tsv(path)
    required = {'conductance', 'module', 'genes_tested', 'n_genes_tested', 'note'}
    required.update(f'{field}_{age}w' for age in AGES for field in (
        'raw_count_proxy', 'proxy_status', 'cpm_primary_factor', 'cpm_status',
        'included_genes', 'missing_raw_counts', 'excluded'))
    if not required.issubset(columns):
        raise ValueError(f'v1 factor file missing required columns: {sorted(required-set(columns))}')
    by_conductance = {}
    for row in rows:
        key = row['conductance']
        if key in by_conductance:
            raise ValueError(f'duplicate v1 factor row: {key}')
        by_conductance[key] = row
    return columns, by_conductance


def complete_row_to_gene(source, conductance, module, gene, tested, source_file=COMPLETE_REL):
    if (source.get('family') != 'whole_transcriptome_secondary'
            or source.get('population') != 'Purkinje cell'
            or source.get('dataset_symbol') != gene
            or not source.get('stable_id')):
        raise ValueError(f'invalid secondary source identity for {gene}: {source}')
    estimate = number(source.get('log2FC_decision'))
    low = number(source.get('CI_lower_decision'))
    high = number(source.get('CI_upper_decision'))
    wt_raw = number(source.get('WT_raw_sum'))
    ki_raw = number(source.get('KI_raw_sum'))
    if wt_raw is not None and wt_raw < 0:
        raise ValueError(f'negative WT raw sum for secondary {gene} at {source["age_weeks"]}w')
    if ki_raw is not None and ki_raw < 0:
        raise ValueError(f'negative KI raw sum for secondary {gene} at {source["age_weeks"]}w')
    failure = source.get('failure', '').strip()
    problem = 'NO_ESTIMATE' if failure or estimate is None else None
    ratio = ratio_low = ratio_high = None
    if problem is None:
        try:
            ratio = 2.0 ** estimate
            ratio_low = 2.0 ** low if low is not None else None
            ratio_high = 2.0 ** high if high is not None else None
        except OverflowError as error:
            raise ValueError(f'ratio overflow for secondary {gene} at {source["age_weeks"]}w') from error
    return {
        'conductance': conductance, 'module': module, 'gene': gene, 'role': 'pore',
        'tested': True, 'age_weeks': age_int(source['age_weeks']),
        'decision_source': source.get('decision_source', ''),
        'n_WT': source.get('n_WT', ''), 'WT_raw_sum': wt_raw, 'KI_raw_sum': ki_raw,
        'log2FC_decision': estimate, 'ratio': ratio, 'ratio_ci_low': ratio_low,
        'ratio_ci_high': ratio_high,
        'raw_count_ratio': ki_raw / wt_raw if wt_raw is not None and wt_raw > 0 and ki_raw is not None else None,
        'problem': problem, 'stable_id': source.get('stable_id', ''),
        'source_file': source_file, 'source_family': 'whole_transcriptome_secondary',
    }


def add_source_labels(gene_rows, release_rows, source_file=RELEASE_REL):
    release_lookup = {}
    for row in release_rows:
        if row.get('population') != 'purkinje':
            continue
        key = (row['dataset_symbol'], age_int(row['age_weeks']))
        if key in release_lookup:
            raise ValueError(f'duplicate release symbol-age: {key}')
        release_lookup[key] = row
    for row in gene_rows:
        source = release_lookup.get((row['gene'], int(row['age_weeks'])))
        row['stable_id'] = source.get('stable_id', '') if source else ''
        row['source_file'] = source_file
        row['source_family'] = source.get('family', 'core_primary') if source else 'core_primary'


def normalized_fields(gene_rows, problems, conductance, module, pore_genes, tested_genes, age):
    tested = sorted(set(pore_genes) & set(tested_genes))
    matches = [r for r in gene_rows if r['conductance'] == conductance
               and r['role'] == 'pore' and r['tested'] is True
               and int(r['age_weeks']) == age and r['gene'] in tested]
    by_gene = {}
    for row in matches:
        if row['gene'] in by_gene:
            raise ValueError(f'duplicate normalized gene row: {conductance}/{row["gene"]}/{age}')
        by_gene[row['gene']] = row
    missing_problem_genes = sorted({p['gene'] for p in problems
                                     if p.get('conductance') == conductance
                                     and int(p.get('age', -1)) == age
                                     and p.get('gene') in tested})
    ratio_rows = []
    excluded = set(missing_problem_genes)
    excluded_weight = set()
    for gene in tested:
        row = by_gene.get(gene)
        ratio = number(row.get('ratio')) if row else None
        problem = row.get('problem') if row else 'MISSING_AGE'
        if row is None or problem in ('NO_ESTIMATE', 'EXCLUDED_UNTESTED_OR_AUX') or ratio is None or ratio <= 0:
            excluded.add(gene)
            continue
        ratio_rows.append(row)
        wt = number(row.get('WT_raw_sum'))
        if wt is None or wt < 0:
            excluded_weight.add(gene)
    weighted_rows = [r for r in ratio_rows
                     if number(r.get('WT_raw_sum')) is not None and number(r.get('WT_raw_sum')) >= 0]
    weight_sum = sum(number(r['WT_raw_sum']) for r in weighted_rows)
    if not tested:
        f_norm = None
        f_status = 'NO_TESTED_GENE'
        geometric = None
        geometric_status = 'NO_TESTED_GENE'
    elif not ratio_rows:
        f_norm = None
        f_status = 'NO_ESTIMATE_AT_AGE'
        geometric = None
        geometric_status = 'NO_ESTIMATE_AT_AGE'
    else:
        geometric = math.exp(sum(math.log(number(r['ratio'])) for r in ratio_rows) / len(ratio_rows))
        geometric_status = 'OK_GEOMETRIC_MEAN_SENSITIVITY'
        if not weighted_rows:
            f_norm = None
            f_status = 'MISSING_WT_WEIGHTS'
        elif weight_sum <= 0:
            f_norm = None
            f_status = 'ZERO_REFERENCE'
        else:
            f_norm = sum(number(r['WT_raw_sum']) * number(r['ratio']) for r in weighted_rows) / weight_sum
            f_status = 'OK_NORMALIZED_EFFECT_WEIGHTED_RAW_WEIGHT_PROXY'
    sources = sorted({r['source_family'] for r in ratio_rows if r.get('source_family')})
    return {
        f'f_norm_{age}w': f_norm,
        f'f_norm_status_{age}w': f_status,
        f'f_norm_weight_sum_{age}w': weight_sum if weighted_rows else None,
        f'f_norm_included_genes_{age}w': ','.join(sorted(r['gene'] for r in weighted_rows)) or '-',
        f'f_norm_excluded_genes_{age}w': ','.join(sorted(excluded | excluded_weight)) or '-',
        f'f_norm_source_families_{age}w': ','.join(sources) or '-',
        f'geometric_mean_sensitivity_{age}w': geometric,
        f'geometric_mean_status_{age}w': geometric_status,
        f'geometric_mean_included_genes_{age}w': ','.join(sorted(r['gene'] for r in ratio_rows)) or '-',
        f'geometric_mean_excluded_genes_{age}w': ','.join(sorted(excluded)) or '-',
    }


def recompute_norm_from_gene_rows(gene_rows, problems, table, tested_by_module):
    recomputed = {}
    for conductance, _, module, pore, _, _ in table.FAMILIES:
        for age in AGES:
            recomputed[(conductance, age)] = normalized_fields(
                gene_rows, problems, conductance, module, pore,
                tested_by_module.get(module, set()), age)
    return recomputed


def same_number(left, right, tolerance=1e-12):
    a, b = number(left), number(right)
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def build_products(release_rows, complete_core, complete_secondary, freeze,
                   v1_factor_columns, v1_factor_file_rows, complete_sha=COMPLETE_SHA,
                   release_source_file=RELEASE_REL, complete_source_file=COMPLETE_REL):
    """Pure derivation from supplied rows; fixtures call this with synthetic-only inputs."""
    table = v1.get_tables()
    tested_by_module = v1.tested_symbols(freeze)
    base_gene_rows, calculated_v1_factors, problems = v1.convert(release_rows, freeze)
    calc_by_conductance = {r['conductance']: r for r in calculated_v1_factors}
    if set(calc_by_conductance) != set(v1_factor_file_rows):
        raise ValueError('calculated v1 conductances differ from the pinned v1 factor file')
    if tuple(v1_factor_columns) == ():
        raise ValueError('v1 factor file has no columns')
    for conductance, computed in calc_by_conductance.items():
        recorded = v1_factor_file_rows[conductance]
        for field in v1_factor_columns:
            if row_text(computed.get(field)) != recorded.get(field, ''):
                raise ValueError(f'v1 factor reproduction mismatch: {conductance}/{field}: '
                                 f'{row_text(computed.get(field))!r} != {recorded.get(field)!r}')

    add_source_labels(base_gene_rows, release_rows, release_source_file)
    output_gene_rows = [dict(row) for row in base_gene_rows]
    problems = [dict(p) for p in problems]

    bca_rows = [family for family in table.FAMILIES if family[0] == 'b_Ca']
    if len(bca_rows) != 1:
        raise ValueError(f'expected exactly one b_Ca family, found {len(bca_rows)}')
    bca_conductance, _, bca_module, bca_pore, _, _ = bca_rows[0]
    tested_bca = sorted(set(bca_pore) & set(tested_by_module.get(bca_module, set())))
    release_symbols = {r['dataset_symbol'] for r in release_rows if r.get('population') == 'purkinje'}
    supplemental_symbols = set(tested_bca) - release_symbols
    if not supplemental_symbols:
        raise ValueError('no tested b_Ca genes are absent from the pinned S08 release')

    for gene in sorted(supplemental_symbols):
        for age in AGES:
            source = complete_secondary.get((gene, 'Purkinje cell', age))
            if source is None:
                continue
            new_row = complete_row_to_gene(source, bca_conductance, bca_module, gene, True,
                                            complete_source_file)
            key = (bca_conductance, gene, age)
            replaced = False
            for index, old in enumerate(output_gene_rows):
                if (old['conductance'], old['gene'], int(old['age_weeks'])) == key:
                    output_gene_rows[index] = new_row
                    replaced = True
                    break
            if not replaced:
                output_gene_rows.append(new_row)
            problems = [p for p in problems if not (
                p.get('conductance') == bca_conductance and p.get('gene') == gene
                and int(p.get('age', -1)) == age)]
            if new_row['problem'] == 'NO_ESTIMATE':
                problems.append({'conductance': bca_conductance, 'gene': gene,
                                 'age': age, 'problem': 'NO_ESTIMATE'})

    output_factor_rows = [dict(v1_factor_file_rows[c]) for c in v1_factor_file_rows]
    factor_by_conductance = {r['conductance']: r for r in output_factor_rows}
    normalized = recompute_norm_from_gene_rows(output_gene_rows, problems, table, tested_by_module)
    for key, fields in normalized.items():
        factor = factor_by_conductance[key[0]]
        for field, value in fields.items():
            factor[field] = value

    # Verify the v1 raw proxy independently from the final per-gene table,
    # using only the release/core rows so b_Ca secondary rows cannot alter it.
    raw_proxy_checks = 0
    for conductance, factor in factor_by_conductance.items():
        for age in AGES:
            at_age = [r for r in output_gene_rows
                      if r['conductance'] == conductance and r['role'] == 'pore'
                      and r['tested'] is True and int(r['age_weeks']) == age
                      and r.get('source_family') == 'core_primary' and r.get('problem') is None]
            wt = [number(r.get('WT_raw_sum')) for r in at_age]
            ki = [number(r.get('KI_raw_sum')) for r in at_age]
            raw_name = f'raw_count_proxy_{age}w'
            expected = factor[raw_name]
            if any(value is None for value in wt + ki):
                if expected != '':
                    raise ValueError(f'v1 raw proxy unexpectedly present with missing sums: {conductance}/{age}')
            elif at_age and sum(wt) > 0:
                observed = sum(ki) / sum(wt)
                if not same_number(observed, expected):
                    raise ValueError(f'raw proxy no longer recomputes: {conductance}/{age}')
            elif expected != '':
                raise ValueError(f'v1 raw proxy unexpectedly present with no reference: {conductance}/{age}')
            raw_proxy_checks += 1

    # No failed or no-estimate row can appear in either normalized aggregate.
    excluded_check = True
    for factor in output_factor_rows:
        for age in AGES:
            included = set(filter(None, factor[f'f_norm_included_genes_{age}w'].split(',')))
            gm_included = set(filter(None, factor[f'geometric_mean_included_genes_{age}w'].split(',')))
            bad = {r['gene'] for r in output_gene_rows if r['conductance'] == factor['conductance']
                   and int(r['age_weeks']) == age
                   and r.get('problem') in ('NO_ESTIMATE', 'EXCLUDED_UNTESTED_OR_AUX')}
            if included & bad or gm_included & bad:
                excluded_check = False
    if not excluded_check:
        raise ValueError('a failed/no-estimate gene entered a normalized family factor')

    buffer_rows = [r for r in output_gene_rows if r['conductance'] == 'b_Ca'
                   and r.get('source_family') == 'whole_transcriptome_secondary']
    aux_rows = [r for r in output_gene_rows if r['role'] == 'auxiliary']
    expected_buffer_rows = len(supplemental_symbols) * len(AGES)
    checks = [
        {'check': 'every release row used matches the complete S08 core row for effects, CIs and raw sums',
         'passed': True, 'n_rows': len(release_rows)},
        {'check': 'normalized factors are recomputable from emitted per-gene rows',
         'passed': True, 'n_factors': len(normalized)},
        {'check': 'v1 raw-count proxy fields reproduce the pinned v1 factor file exactly',
         'passed': True, 'n_factor_age_cells': raw_proxy_checks},
        {'check': 'failed/no-estimate genes are excluded from normalized aggregates and listed',
         'passed': excluded_check},
        {'check': 'secondary b_Ca rows carry their source file and family',
         'passed': all(r.get('source_file') == complete_source_file and
                       r.get('source_family') == 'whole_transcriptome_secondary'
                       for r in buffer_rows) and len(buffer_rows) == expected_buffer_rows,
         'n_rows': len(buffer_rows), 'expected_rows': expected_buffer_rows},
        {'check': 'exact CPM primary factor remains unavailable from the aggregate release',
         'passed': all(factor_by_conductance[c].get(f'cpm_primary_factor_{a}w', '') == ''
                       for c in factor_by_conductance for a in AGES)},
    ]
    if any(not c['passed'] for c in checks):
        raise ValueError(f'B3 v02 self-check failed: {checks}')
    problems.sort(key=lambda p: (p.get('conductance', ''), p.get('gene', ''),
                                 int(p.get('age', 0)), p.get('problem', '')))
    return {
        'gene_rows': output_gene_rows,
        'factor_rows': output_factor_rows,
        'factor_columns': tuple(v1_factor_columns) + tuple(
            f'{suffix}_{age}w' for age in AGES for suffix in NORM_SUFFIXES),
        'auxiliary_rows': aux_rows,
        'buffer_rows': buffer_rows,
        'problems': problems,
        'checks': checks,
        'tested_by_module': tested_by_module,
        'table': table,
        'complete_sha256': complete_sha,
        'release_source_file': release_source_file,
        'complete_source_file': complete_source_file,
    }


def tsv_rows_to_dicts(path):
    _, rows = read_tsv(path)
    for row in rows:
        for key in ('tested',):
            if row.get(key) in ('True', 'False'):
                row[key] = row[key] == 'True'
        for key in ('age_weeks',):
            if row.get(key):
                row[key] = int(row[key])
    return rows


def recompute_serialized_factors(gene_path, factor_path, problems, table, tested_by_module):
    genes = tsv_rows_to_dicts(gene_path)
    factors = tsv_rows_to_dicts(factor_path)
    calculated = recompute_norm_from_gene_rows(genes, problems, table, tested_by_module)
    for factor in factors:
        for age in AGES:
            expected = calculated[(factor['conductance'], age)]
            for suffix in ('f_norm', 'geometric_mean_sensitivity'):
                field = f'{suffix}_{age}w'
                if not same_number(factor.get(field), expected.get(field)):
                    raise ValueError(f'serialized gene rows do not reproduce {field} for {factor["conductance"]}')
            for suffix in ('f_norm_status', 'f_norm_included_genes', 'f_norm_excluded_genes',
                           'f_norm_source_families', 'geometric_mean_status',
                           'geometric_mean_included_genes', 'geometric_mean_excluded_genes'):
                field = f'{suffix}_{age}w'
                if factor.get(field, '') != row_text(expected.get(field)):
                    raise ValueError(f'serialized gene rows do not reproduce {field} for {factor["conductance"]}')
    return len(calculated)


def build_result(products, input_hashes, serialized_norm_cells):
    checks = [dict(c) for c in products['checks']]
    checks[1]['passed'] = serialized_norm_cells == len(products['factor_rows']) * len(AGES)
    checks[1]['n_factors'] = serialized_norm_cells
    return {
        'version': 'B3_v02_build',
        'mode': 'release_level_derived_factors; no model run',
        'rule': ('F_norm is a descriptive WT_raw_sum-weighted mean of normalized per-gene 2^log2FC_decision effects; '
                 'because WT_raw_sum is unnormalized, F_norm is a raw-weight proxy sensitivity, not a CPM family factor '
                 'or worker scale setting. The unweighted geometric mean is S1 sensitivity. The v1 raw-count proxy '
                 'remains unchanged and CPM_PRIMARY_FACTOR_UNAVAILABLE_FROM_RELEASE remains explicit. '
                 'b_Ca secondary estimates are sensitivity-only.'),
        'input_sha256': input_hashes,
        'complete_s08_sha256': products['complete_sha256'],
        'n_gene_rows': len(products['gene_rows']),
        'n_factors': len(products['factor_rows']),
        'n_buffer_sensitivity_rows': len(products['buffer_rows']),
        'n_problems': len(products['problems']),
        'problems': products['problems'],
        'n_checks': len(checks), 'n_pass': sum(bool(c['passed']) for c in checks),
        'all_pass': all(bool(c['passed']) for c in checks), 'checks': checks,
        'model_loaded': False, 'simulation_run': False,
    }


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 5:
        raise SystemExit('expected: <release.tsv.gz> <FREEZE.json> <complete.tsv.gz> <v1 CONDUCTANCE_FACTORS.tsv> <new output dir>')
    release_path, freeze_path, complete_path, v1_factor_path, out = map(Path, args)
    if out.exists():
        raise FileExistsError(f'append-only: output path already exists: {out}')

    # Check converter/table identity before reading any expression or FREEZE input.
    if sha256_file(SCRIPTS / 'p4_t2_scale_factors_v01.py') != V1_SCRIPT_SHA:
        raise ValueError('v1 converter dependency differs from its pin')
    if sha256_file(v1.TABLES_PATH) != TABLES_SHA:
        raise ValueError('family table module differs from its pin')
    paths = {'release': release_path, 'freeze': freeze_path,
             'complete_results': complete_path, 'v1_factors': v1_factor_path}
    expected = {'release': RELEASE_SHA, 'freeze': FREEZE_SHA,
                'complete_results': COMPLETE_SHA, 'v1_factors': V1_FACTORS_SHA}
    input_hashes = {}
    for name, path in paths.items():
        got = sha256_file(path)
        input_hashes[name] = got
        if got != expected[name]:
            raise ValueError(f'{name} differs from its pin: {got}')

    # The table loader verifies the module bytes before parsing FREEZE/release.
    table, table_sha = v1.load_tables({'tables': TABLES_SHA})
    v1.tables = table
    release_rows = v1.read_release(release_path)
    freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
    tested = v1.tested_symbols(freeze)
    bca = [f for f in table.FAMILIES if f[0] == 'b_Ca']
    if len(bca) != 1:
        raise ValueError(f'expected exactly one b_Ca family, found {len(bca)}')
    _, _, bca_module, bca_pore, _, _ = bca[0]
    released_symbols = {r['dataset_symbol'] for r in release_rows if r['population'] == 'purkinje'}
    supplemental = set(bca_pore) & set(tested.get(bca_module, set())) - released_symbols
    complete_core, complete_secondary = read_complete_index(complete_path, release_rows, supplemental)
    rows_verified = check_release_complete_identity(release_rows, complete_core)
    v1_factor_columns, v1_factor_rows = read_v1_factor_file(v1_factor_path)
    products = build_products(release_rows, complete_core, complete_secondary, freeze,
                              v1_factor_columns, v1_factor_rows, input_hashes['complete_results'])
    products['checks'][0]['n_rows'] = rows_verified
    products['checks'][4]['n_rows'] = len(products['buffer_rows'])

    # Exclusive-create the new v02 destination only after all pins and source checks pass.
    out.mkdir(parents=True, exist_ok=False)
    gene_path = out / 'GENE_RATIOS.tsv'
    factor_path = out / 'CONDUCTANCE_FACTORS.tsv'
    aux_path = out / 'AUXILIARY.tsv'
    buffer_path = out / 'BUFFER_SENSITIVITY.tsv'
    write_tsv_exclusive(gene_path, products['gene_rows'], GENE_COLUMNS)
    write_tsv_exclusive(aux_path, products['auxiliary_rows'], GENE_COLUMNS)
    write_tsv_exclusive(buffer_path, products['buffer_rows'], BUFFER_COLUMNS)
    write_tsv_exclusive(factor_path, products['factor_rows'], products['factor_columns'])
    norm_cells = recompute_serialized_factors(gene_path, factor_path, products['problems'],
                                              products['table'], products['tested_by_module'])
    input_hashes['tables'] = table_sha
    input_hashes['v1_converter'] = V1_SCRIPT_SHA
    result = build_result(products, input_hashes, norm_cells)
    with (out / 'RESULT.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(f"{result['n_pass']} / {result['n_checks']} all_pass {result['all_pass']} "
          f"| factors {result['n_factors']} | b_Ca rows {result['n_buffer_sensitivity_rows']} "
          f"| problems {result['n_problems']}")
    return 0 if result['all_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
