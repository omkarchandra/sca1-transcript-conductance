"""Plot descriptive 18 W KI/WT rate ratios from a population scorer set.

E4 (pre-specified): copy of the Figure 3B print script (a11e4e74) with the
y-axis, x labels (gene/channel) and legend (one entry) changed and the three inputs hash-pinned.
Data, order, jitter and colours unchanged.

python -I -B 06_figures/fig5b_values.py --scored DIR --out NEW_DIR
The real-data run requires a separate review block clearing this script's bytes.
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
from matplotlib.patches import Patch


SINGLES = ('Kcnc3', 'Kcnma1', 'Kcnd3', 'Cacna1a', 'Cacna1g')
GROUPS = ('BK_CaP', 'Kv3_Kv4_CaT', 'BK_CaT', 'Kv3_Kv4_CaP')
CONDITIONS = ('KI_18w_mid',) + tuple(f'KI_18w_mid_restore_{gene}' for gene in SINGLES) + tuple(
    f'KI_18w_mid_grp_{group}' for group in GROUPS)
CACNA1A = frozenset(('KI_18w_mid_restore_Cacna1a', 'KI_18w_mid_grp_BK_CaP',
                     'KI_18w_mid_grp_Kv3_Kv4_CaP'))
PINS = {'SUMMARY.json': '5c0b3a335e960d5c8d772df594406b57f837e6e2e2afa79dab0e507c752311ea',
        'PER_MEMBER.tsv': '148890d00ad1b45d2e5a04842d1083898c56242eedac1817b953ec3ba140d7fa',
        'PER_RUN.tsv': 'b8b9f2a8d3f4e6e6a818626dab06f4bbbb6c2d0910f022b83fae1d7adf57ee5a'}
LABELS = {'KI_18w_mid': 'unrestored KI',
          'KI_18w_mid_restore_Kcnc3': 'restore Kcnc3/Kv3',
          'KI_18w_mid_restore_Kcnma1': 'restore Kcnma1/BK',
          'KI_18w_mid_restore_Kcnd3': 'restore Kcnd3/Kv4',
          'KI_18w_mid_restore_Cacna1a': 'restore Cacna1a/CaP',
          'KI_18w_mid_restore_Cacna1g': 'restore Cacna1g/CaT',
          'KI_18w_mid_grp_BK_CaP': 'Kcnma1/BK + Cacna1a/CaP',
          'KI_18w_mid_grp_Kv3_Kv4_CaT': 'Kcnc3/Kv3 + Kcnd3/Kv4 + Cacna1g/CaT',
          'KI_18w_mid_grp_BK_CaT': 'Kcnma1/BK + Cacna1g/CaT',
          'KI_18w_mid_grp_Kv3_Kv4_CaP': 'Kcnc3/Kv3 + Kcnd3/Kv4 + Cacna1a/CaP'}
Y_LABEL = 'KI-scaled / WT settings, firing rate'
LEGEND = 'orange = Cacna1a/CaP put back to WT in that run'
FIGURE = 'firing_rate_ratio_restorations.png'
SVG_FIGURE = 'firing_rate_ratio_restorations.svg'
SUBTITLE = ('Descriptive, not pre-registered. Unrestored KI already met the reversal rule; '
            'this panel shows the rate change only.')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def table(path, columns):
    with path.open('r', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not set(columns).issubset(reader.fieldnames or ()):
            raise ValueError(f'{path}: missing columns {set(columns) - set(reader.fieldnames or ())}')
        return list(reader)


def read_numbers(scored):
    paths = {name: scored / name for name in ('SUMMARY.json', 'PER_MEMBER.tsv', 'PER_RUN.tsv')}
    for name, path in paths.items():  # E4: inputs pinned to the 840-run scoring of record
        if digest(path) != PINS[name]:
            raise ValueError(f'{name}: input hash differs from the pinned scoring')
    summary = json.loads(paths['SUMMARY.json'].read_text(encoding='utf-8'))
    members = table(paths['PER_MEMBER.tsv'], ('member_id', 'kept', 'condition', 'rate_ratio'))
    table(paths['PER_RUN.tsv'], ('member_id', 'condition', 'rate_hz'))
    if summary.get('n_kept') != 14:
        raise ValueError('expected exactly 14 kept members in scorer summary')
    kept_ids = {int(row['member_id']) for row in members if row['kept'] == 'True'}
    if len(kept_ids) != 14:
        raise ValueError('expected exactly 14 distinct kept member IDs')
    plotted = []
    for condition in CONDITIONS:
        rows = [row for row in members if row['condition'] == condition and row['kept'] == 'True']
        ids = [int(row['member_id']) for row in rows]
        if len(ids) != 14 or set(ids) != kept_ids or len(ids) != len(set(ids)):
            raise ValueError(f'{condition}: expected one comparison for each of 14 kept members')
        points = []
        for row in sorted(rows, key=lambda r: int(r['member_id'])):
            value = float(row['rate_ratio'])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'{condition}: invalid rate ratio')
            points.append(dict(member_id=int(row['member_id']), rate_ratio=value))
        recomputed = median(point['rate_ratio'] for point in points)
        if condition == 'KI_18w_mid':
            reported = summary['predictions'][condition]['median_rate_ratio']
        elif '_restore_' in condition:
            reported = summary['reversal'][condition]['median_rate_ratio']
        else:
            reported = summary['grp_reversal'][condition]['median_rate_ratio']
        if not math.isclose(recomputed, float(reported), rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError(f'{condition}: scorer summary median differs from PER_MEMBER')
        plotted.append(dict(condition=condition, includes_Cacna1a=condition in CACNA1A,
                            n=len(points), median_rate_ratio=recomputed, points=points))
    return plotted, {name: digest(path) for name, path in paths.items()}


def draw(path, plotted):
    fig, ax = plt.subplots(figsize=(3.46, 4.8), constrained_layout=False)
    dark, light = '#163b66', '#f4a261'
    ax.axhspan(.75, 1.25, color='#ced4da', alpha=.35, label='reversal-rule rate tolerance')
    ax.axhline(1, color='#495057', linewidth=.6)
    labels = []
    for i, entry in enumerate(plotted):
        color = light if entry['includes_Cacna1a'] else dark
        for rank, point in enumerate(entry['points']):
            ax.scatter(i + (rank - 6.5) * .012, point['rate_ratio'], s=6,
                       facecolor=color, edgecolor='#1f2933', linewidth=.2, zorder=3)
        assert entry['median_rate_ratio'] == median(p['rate_ratio'] for p in entry['points'])
        ax.plot([i - .2, i + .2], [entry['median_rate_ratio']] * 2,
                color='#212529', linewidth=1.2, zorder=4)
        labels.append(LABELS[entry['condition']])
    ax.set_xticks(range(len(plotted)), labels, fontsize=6, rotation=55, ha='right', rotation_mode='anchor')
    ax.set_xlim(-.55, len(plotted) - .45)
    ax.set_ylabel(Y_LABEL, fontsize=7)
    ax.tick_params(axis='y', labelsize=6)
    fig.suptitle('KI firing-rate change and restorations (18 W)', fontsize=8, y=.98)
    fig.legend(handles=[Patch(color=light, label=LEGEND)],
               loc='lower center', bbox_to_anchor=(.55, .005), ncol=1, frameon=False, fontsize=6)
    fig.subplots_adjust(left=.2, right=.98, bottom=.36, top=.91)
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
                  title='KI firing-rate change and restorations (18 W)', subtitle=SUBTITLE, y_label=Y_LABEL,
                  x_labels=LABELS, legend=[LEGEND],
                  source_script_sha256='a11e4e74d6fba84084b809b8b7ec82e574bf1a17d8f16375b7b8dd08fb4232ad',
                  band=dict(lower=.75, upper=1.25, reference=1.0,
                            label='reversal-rule rate tolerance'), numbers_plotted=plotted)
    with (args.out / 'FIGURE_RECORD.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(record, indent=2, allow_nan=False) + '\n')
    with (args.out / 'SHA256SUMS').open('x', encoding='utf-8') as stream:
        for name in (FIGURE, SVG_FIGURE, 'FIGURE_RECORD.json'):
            stream.write(f'{digest(args.out / name)}  {name}\n')


if __name__ == '__main__':
    main()

#%% debug
