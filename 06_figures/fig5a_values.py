"""Descriptive per-member KI minus own-WT ISI-CV panel from scorer outputs.

E4 (pre-specified): copy of the Figure 3A print script (c6571a6f) with the
y-axis relabelled and the three inputs hash-pinned. Data, positions, jitter and styling unchanged.

python -I -B 06_figures/fig5a_values.py
    --scored DIR --out NEW_DIR
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['svg.hashsalt'] = 'phase4'
matplotlib.use('Agg')
import matplotlib.pyplot as plt


CONDITIONS = ('KI_5w_mid',) + tuple(
    f'KI_{age}w_{level}' for age in (18, 24, 30) for level in ('low', 'mid', 'high'))
POSITIONS = (0, 1.5, 2.5, 3.5, 5, 6, 7, 8.5, 9.5, 10.5)
AGES = (5, 18, 24, 30)
PINS = {'SUMMARY.json': '5c0b3a335e960d5c8d772df594406b57f837e6e2e2afa79dab0e507c752311ea',
        'PER_MEMBER.tsv': '148890d00ad1b45d2e5a04842d1083898c56242eedac1817b953ec3ba140d7fa',
        'PER_RUN.tsv': 'b8b9f2a8d3f4e6e6a818626dab06f4bbbb6c2d0910f022b83fae1d7adf57ee5a'}
Y_LABEL = 'Δ ISI CV (KI-scaled − WT settings)'
FIGURE = 'fig2_delta_cv_by_condition.png'
SVG_FIGURE = 'fig2_delta_cv_by_condition.svg'
TITLE = 'KI change in firing irregularity vs own WT'
SUBTITLE = ('Each dot is one model cell; grey ticks mark the increase it needed to count as more irregular. '
            'Pre-registered prediction 1: more irregular at 18–30 W.')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path, required):
    with path.open('r', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or ()):
            raise ValueError(f'{path}: missing columns {set(required) - set(reader.fieldnames or ())}')
        return list(reader)


def boolean(text):
    if text not in ('True', 'False'):
        raise ValueError(f'invalid boolean: {text!r}')
    return text == 'True'


def finite(text):
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f'non-finite number: {text!r}')
    return value


def read_numbers(scored):
    paths = {name: scored / name for name in ('SUMMARY.json', 'PER_MEMBER.tsv', 'PER_RUN.tsv')}
    for name, path in paths.items():  # E4: inputs pinned to the 840-run scoring of record
        if digest(path) != PINS[name]:
            raise ValueError(f'{name}: input hash differs from the pinned scoring')
    summary = json.loads(paths['SUMMARY.json'].read_text(encoding='utf-8'))
    members = table(paths['PER_MEMBER.tsv'],
                    ('member_id', 'kept', 'condition', 'cv_increase', 'cv_margin', 'cv_increased'))
    table(paths['PER_RUN.tsv'], ('member_id', 'condition', 'isi_cv'))
    if summary['n_kept'] != 14:
        raise ValueError('expected 14 kept members in scorer summary')
    plotted = []
    reference_ids = None
    for condition, x in zip(CONDITIONS, POSITIONS):
        rows = [r for r in members if r['condition'] == condition and boolean(r['kept'])]
        ids = [int(r['member_id']) for r in rows]
        if len(rows) != 14 or len(set(ids)) != 14:
            raise ValueError(f'{condition}: expected 14 unique kept members')
        if reference_ids is None:
            reference_ids = set(ids)
        if set(ids) != reference_ids:
            raise ValueError(f'{condition}: kept member IDs differ')
        points = [dict(member_id=int(r['member_id']), cv_increase=finite(r['cv_increase']),
                       cv_margin=finite(r['cv_margin']), cv_increased=boolean(r['cv_increased']))
                  for r in sorted(rows, key=lambda r: int(r['member_id']))]
        if any(point['cv_margin'] < 0 for point in points):
            raise ValueError(f'{condition}: negative threshold')
        recomputed = median(point['cv_increase'] for point in points)
        n_increased = sum(point['cv_increased'] for point in points)
        section = 'predictions' if condition.endswith('_mid') else 'bounds'
        entry = summary[section][condition]
        if not math.isclose(recomputed, float(entry['median_cv_increase']), rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError(f'{condition}: median differs from PER_MEMBER')
        if n_increased != entry['n_cv_increased'] or len(points) != entry['n_kept_members_with_run']:
            raise ValueError(f'{condition}: count differs from SUMMARY')
        plotted.append(dict(condition=condition, age=int(condition.split('_')[1][:-1]),
                            level=condition.rsplit('_', 1)[-1], x=x, n=len(points),
                            n_increased=n_increased, median_cv_increase=recomputed, points=points))
    return plotted, {name: digest(path) for name, path in paths.items()}


def draw(path, plotted):
    fig, ax = plt.subplots(figsize=(3.46, 4.6))
    palette = {5: '#174e8c', 18: '#087f8c', 24: '#7b4ea3', 30: '#c35b20'}
    for age, lo, hi in ((5, -.46, .46), (18, 1.04, 3.96),
                        (24, 4.54, 7.46), (30, 8.04, 10.96)):
        if age in (18, 30):
            ax.axvspan(lo, hi, color='#edf2f5', zorder=0)
    ax.axhline(0, color='#343a40', linewidth=.55, zorder=1)
    for entry in plotted:
        x = entry['x']
        for rank, point in enumerate(entry['points']):
            jitter = (rank - 6.5) * .025
            px = x + jitter
            ax.plot([px - .012, px + .012], [point['cv_margin']] * 2,
                    color='#777e85', linewidth=.55, zorder=2)
            ax.scatter(px, point['cv_increase'], s=6, color=palette[entry['age']],
                       edgecolors='#22303a', linewidths=.2, zorder=3)
        assert entry['median_cv_increase'] == median(p['cv_increase'] for p in entry['points'])
        ax.plot([x - .24, x + .24], [entry['median_cv_increase']] * 2,
                color='#1c2530', linewidth=1.25, zorder=4)
        ax.text(x, 1.01, f"{entry['n_increased']}/{entry['n']}",
                ha='center', va='bottom', fontsize=6, transform=ax.get_xaxis_transform())
    ax.set_xticks(POSITIONS, ['mid'] + ['low', 'mid', 'high'] * 3, fontsize=6)
    for age, x in ((5, 0), (18, 2.5), (24, 6), (30, 9.5)):
        ax.text(x, -.15, f'{age} W', ha='center', va='top', fontsize=6,
                transform=ax.get_xaxis_transform())
    ax.set_xlim(-.55, 11.05)
    ax.set_ylabel(Y_LABEL, fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    fig.suptitle(TITLE, fontsize=8, y=.98)
    fig.subplots_adjust(left=.2, right=.98, top=.92, bottom=.22)
    fig.savefig(path.with_name(SVG_FIGURE), metadata={'Date': None})
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scored', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    plotted, hashes = read_numbers(args.scored)
    args.out.mkdir(parents=False, exist_ok=False)
    draw(args.out / FIGURE, plotted)
    record = dict(script_sha256=digest(__file__), input_sha256=hashes,
                  title=TITLE, subtitle=SUBTITLE, y_label=Y_LABEL, source_script_sha256='c6571a6ff0b3ca693e8d244c8c39d2a258949982374d1fab7bd6d85bd0624294',
                  reference_line=0, x_conditions=list(CONDITIONS),
                  numbers_plotted=plotted)
    with (args.out / 'FIGURE_RECORD.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(record, indent=2, allow_nan=False) + '\n')
    with (args.out / 'SHA256SUMS').open('x', encoding='utf-8') as stream:
        for name in (FIGURE, SVG_FIGURE, 'FIGURE_RECORD.json'):
            stream.write(f'{digest(args.out / name)}  {name}\n')


if __name__ == '__main__':
    main()

#%% debug
