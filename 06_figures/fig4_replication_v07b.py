"""Plot the nine-gene replication estimates and 95 % CIs as Figure 4, as drawn for the manuscript."""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

# The bars, CIs, dataset titles and counts come from 06_figures/replication_values.py (pinned below); the table written
# with the figure must equal results/replication/PLOTTED_BARS.tsv byte for byte. Labels added after the draw:
#   - no top-line title and no footnote; shared y-axis label "log2 fold change";
#   - x labels: channel or protein name, then the gene in brackets (two lines, 8 pt); group headers 'Five channel
#     transcripts' / 'Four other transcripts'; each row's count of negative estimates (5 channels, 9 transcripts);
#   - row titles of the 12-week held-out row and the ATXN1[82Q] row as in the manuscript (DISPLAY);
#   - one asterisk per transcript with q < 0.05: held-out row from q_core of data/snrnaseq/E3_LIMMA_RESULTS.tsv; the
#     ATXN1[82Q], SCA7 and human rows from the stored adj.P.Val of the a1_limma.R runs (NINE_GENES.tsv; topTable over all
#     genes, BH); the bulk row is not marked (authors' values). The marked sets are checked before anything is drawn;
#   - a two-line note under the x labels. No artist of the helper is moved and the canvas is unchanged.
# The earlier 06_figures/fig4_replication.py draws the same bars with the earlier labels.
# Run: python -I -B 06_figures/fig4_replication_v07b.py
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(_release_paths.expand('${REPO_ROOT}'))
HERE = Path(__file__).resolve().parent
SOURCE = ROOT / '06_figures/replication_values.py'
SOURCE_SHA = '68d63352f564fba72c53431733eae7c1d2d8e74c4e262639fa97b77f850f33a5'
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/fig4_replication_v07b'))
NAMES = {'Kcnma1': 'BK', 'Kcnc3': 'Kv3.3', 'Kcnd3': 'Kv4.3', 'Cacna1a': 'Cav2.1', 'Cacna1g': 'Cav3.1',
         'Itpr1': 'IP3R1', 'Atp2a3': 'SERCA3', 'Kcnmb4': 'BK β4', 'Cacna2d2': 'α2δ2'}
HEADERS = {'Five channel genes': 'Five channel transcripts', 'Four other genes': 'Four other transcripts'}
TICK_PT, MIN_GAP_CM = 8, 0.1
DISPLAY = {'Discovery Purkinje nuclei, 12 weeks (held-out mice)': 'Knock-in Purkinje nuclei, 12 weeks (held-out mice)',
           'ATXN1[82Q] Purkinje nuclei, 12 weeks': 'ATXN1[82Q] Purkinje nuclei, 12 weeks (no single transcript q < 0.05)'}
Q_SOURCE = ROOT / 'data/snrnaseq/E3_LIMMA_RESULTS.tsv'
Q_EXPECTED = {'Kcnma1': 0.0458, 'Kcnc3': 0.0006, 'Kcnd3': 0.0244, 'Cacna1a': 0.0024, 'Cacna1g': 0.0197,
              'Itpr1': 0.1048, 'Atp2a3': 0.0066, 'Kcnmb4': 0.0167, 'Cacna2d2': 0.0072}
NOTE = ("* q < 0.05 (Benjamini–Hochberg); the bulk row is not marked (authors' values).\n"
        'Whiskers: 95 % confidence intervals, not corrected for the number of genes tested.')
NOTE_TWO_LINES = NOTE.replace('Hochberg); ', 'Hochberg);\n')
A1_LIMMA_AS_RUN = '10c810956e27ca042bd8c927e01f9b6f559b3d323dec270f068a91fb4d77f0c9'
EXPECTED_MARKED = {
    'ATXN1[82Q] Purkinje nuclei, 12 weeks': set(),
    'Human SCA1 Purkinje nuclei, 13 donors (primary)': set(),
    'SCA7 Purkinje nuclei, 5 weeks': {'Kcnma1', 'Cacna1a', 'Itpr1', 'Kcnmb4'},
    'SCA7 Purkinje nuclei, 8.5 weeks': {'Kcnma1', 'Kcnc3', 'Cacna1a', 'Cacna1g', 'Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2'}}
MIN_CLEAR_PX300, MIN_NOTE_GAP_CM, NOTE_GAP_CM = 2, 0.1, 0.15


def main():
    assert not OUT.exists(), 'exclusive-create output; do not overwrite'
    spec = importlib.util.spec_from_file_location('accepted_replication', SOURCE)
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    assert b.sha(SOURCE) == SOURCE_SHA
    table = b.load()
    accepted_path = ROOT / 'results/replication/PLOTTED_BARS.tsv'
    accepted = b.pd.read_csv(accepted_path, sep='\t')
    assert list(table.gene) == list(accepted.gene) and list(table.dataset) == list(accepted.dataset)
    for col in ('log2FC', 'ci_low', 'ci_high'):
        assert np.allclose(table[col], accepted[col], rtol=0, atol=1e-14)
    # v06: q values of the nine transcripts at 12 weeks, Purkinje cells; stop before drawing on any mismatch.
    q_sha = b.sha(Q_SOURCE)
    lr = [l.rstrip('\n').split('\t') for l in Q_SOURCE.open(encoding='utf-8')]
    lr = [dict(zip(lr[0], r)) for r in lr[1:]]
    q = {}
    for g in b.GENES:
        m = [r for r in lr if r['dataset_symbol'] == g and r['population'] == 'Purkinje cell' and r['age_weeks'] == '12']
        assert len(m) == 1, (g, len(m))
        q[g] = float(m[0]['q_core'])
    assert {g: round(q[g], 4) for g in b.GENES} == Q_EXPECTED, ('q mismatch', q)
    marked = [g for g in b.GENES if q[g] < 0.05]
    assert len(marked) == 8 and 'Itpr1' not in marked, marked
    held = b.DATASETS[0][0]
    n_core = sum(1 for r in lr if r['population'] == 'Purkinje cell' and r['age_weeks'] == '12' and r['q_core'] not in ('', 'NA'))
    # v07: the other single-nucleus rows. q = the stored adj.P.Val of the a1_limma.R run that drew the row (BH over all
    # genes tested); log2FC must equal the plotted bar; the marked set must be the reported one. Stop before drawing otherwise.
    marks, qrows = {held: list(marked)}, []
    for g in b.GENES:
        r0 = accepted[(accepted.dataset == held) & (accepted.gene == g)].iloc[0]
        qrows.append(dict(dataset=held, gene=g, log2FC=float(r0.log2FC), q=q[g], marked='yes' if g in marked else 'no',
                          n_genes_tested='%d (q_core: BH over the core-panel genes tested in Purkinje cells at 12 weeks; e3_limma_12w.R line 146: q_core[core_hit] <- p.adjust(p[core_hit], "BH"))' % n_core, source=str(Q_SOURCE), source_sha256=q_sha))
    q_files = {}
    for title, verdict, source, dataset, age in b.DATASETS:
        if source == 'p2':
            for g in b.GENES:
                r0 = accepted[(accepted.dataset == title) & (accepted.gene == g)].iloc[0]
                qrows.append(dict(dataset=title, gene=g, log2FC=float(r0.log2FC), q='NA', marked='no (authors\' values)',
                                  n_genes_tested='NA', source='authors\' values', source_sha256='NA'))
            continue
        if source not in ('a1', 'expl'):
            continue
        nine_path, nine_sha = b.SOURCES[source]
        read_path, read_sha = b.SOURCES['expl_reading' if source == 'expl' else 'a1_reading']
        assert b.sha(nine_path) == nine_sha and b.sha(read_path) == read_sha, (title, 'source changed')
        q_files[source] = dict(nine_genes=dict(path=str(nine_path), sha256=nine_sha), reading=dict(path=str(read_path), sha256=read_sha))
        t = b.pd.read_csv(nine_path, sep='\t', comment='#')
        sel = t[(t.dataset == dataset) & (t.age == age)]
        sel = sel[sel.gene.isin(b.GENES)].set_index('gene').loc[b.GENES]
        rd = b.pd.read_csv(read_path, sep='\t', comment='#')
        rd = rd[(rd.dataset == dataset) & (rd.age == age)]
        assert len(sel) == 9 and len(rd) == 1, (title, len(sel), len(rd))
        n_tested = int(rd.iloc[0].n_genes_tested)
        assert n_tested > 1000, (title, n_tested)
        acc = accepted[accepted.dataset == title].set_index('gene').loc[b.GENES]
        assert np.allclose(sel.logFC.to_numpy(float), acc.log2FC.to_numpy(float), rtol=0, atol=1e-12), (title, 'log2FC differs from PLOTTED_BARS')
        qv = {g: float(sel.loc[g, 'adj.P.Val']) for g in b.GENES}
        row_marked = [g for g in b.GENES if qv[g] < 0.05]
        assert set(row_marked) == EXPECTED_MARKED[title], (title, row_marked)
        marks[title] = row_marked
        for g in b.GENES:
            qrows.append(dict(dataset=title, gene=g, log2FC=float(sel.loc[g, 'logFC']), q=qv[g], marked='yes' if g in row_marked else 'no',
                              n_genes_tested=str(n_tested), source=str(nine_path), source_sha256=nine_sha))
    assert set(marks) == {held} | set(EXPECTED_MARKED), sorted(marks)
    fig = plt.figure(figsize=(b.WIDTH / 2.54, b.HEIGHT / 2.54))
    b.draw(fig, table)
    # v03: change only right-row labels; accepted bars and CIs untouched.
    row_labels = {'Discovery Purkinje nuclei, 12 weeks (held-out mice)': '5/5, 9/9', 'ATXN1[82Q] Purkinje nuclei, 12 weeks': '5/5, 9/9', 'SCA7 Purkinje nuclei, 8.5 weeks': '5/5, 9/9', 'SCA7 Purkinje nuclei, 5 weeks': '5/5, 9/9', 'Knock-in bulk cerebellum, 12 weeks': '3/5, 6/9', 'Human SCA1 Purkinje nuclei, 13 donors (primary)': '2/5, 5/9'}
    for title, verdict, *_ in b.DATASETS:
        plotted = table[table.dataset == title]
        channels = plotted[plotted.gene.isin(['Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g'])]
        counts = (int((channels.log2FC < 0).sum()), int((plotted.log2FC < 0).sum()))
        assert len(channels) == 5 and len(plotted) == 9, (title, len(channels), len(plotted))
        assert row_labels[title] == f'{counts[0]}/5, {counts[1]}/9', (title, counts, row_labels[title])
        name = next(t for t in fig.texts if t.get_text() == title)
        right = next(t for t in fig.texts if t.get_text() == verdict
                     and t.get_position()[1] == name.get_position()[1])
        right.set_text(row_labels[title])
    # v05: change only the two drawn row titles, after the rows were found by their original titles.
    for old, new in DISPLAY.items():
        hits = [t for t in fig.texts if t.get_text() == old]
        assert len(hits) == 1, (old, len(hits))
        hits[0].set_text(new)
    # v02 label-only surgery after the accepted draw(); no data or geometry artist is touched.
    for text in list(fig.texts):
        # v01: only the panel letter is removed; all data geometry stays as accepted.
        if text.get_text() == 'B':
            text.remove()
        elif text.get_text().startswith('Replication: KI'):
            text.remove()  # top-line title; its quantity moves to the shared y-axis label
        elif text.get_text().startswith('Bars: accepted estimates'):
            text.remove()  # in-image explanatory footnote (colour legends are kept)
    # v04: two-line x labels on the bottom strip (the only one with visible x labels); group headers renamed.
    assert len(fig.axes) == 6
    bottom_ax = fig.axes[5]
    assert [t.get_text() for t in bottom_ax.get_xticklabels()] == list(b.GENES)
    bottom_ax.set_xticks(bottom_ax.get_xticks(), [f'{NAMES[g]}\n({g})' for g in b.GENES], fontsize=TICK_PT)
    for text in fig.texts:
        if text.get_text() in HEADERS:
            text.set_text(HEADERS[text.get_text()])
    # Shared y-axis label, rotated, centred on the six-strip stack (axes span 15.05 to 2.30 cm of 17.0).
    ylabel = fig.text(.45 / b.WIDTH, 8.675 / b.HEIGHT, 'log2 fold change',
                      fontsize=9, rotation=90, ha='center', va='center')
    bulk = fig.axes[4]
    marker = bulk.text(1, .35, '†', fontsize=9, ha='center', va='bottom', color='#333333', zorder=5)
    assert marker.get_visible() and marker.get_text() == '†'
    # v07: one asterisk per marked transcript in its own row, same font size as the bulk-row dagger: at the dagger height
    # for a bar below zero, just above the whisker for a bar above zero.
    row1 = fig.axes[0]
    assert all(ax.get_ylim() == bulk.get_ylim() and ax.get_xlim() == bulk.get_xlim() for ax in fig.axes)
    titles = [d[0] for d in b.DATASETS]
    stars = {}
    for title, gs in marks.items():
        ax = fig.axes[titles.index(title)]
        rows = table[table.dataset == title].set_index('gene')
        for g in gs:
            j, v, hi = b.GENES.index(g), float(rows.loc[g, 'log2FC']), float(rows.loc[g, 'ci_high'])
            if v < 0:
                stars[(title, g)] = ax.text(j, marker.get_position()[1], '*', fontsize=marker.get_fontsize(), ha='center', va='bottom', color='#333333', zorder=5)
            else:
                stars[(title, g)] = ax.annotate('*', xy=(j, hi), xytext=(0, 2), textcoords='offset points', fontsize=marker.get_fontsize(),
                                                ha='center', va='bottom', color='#333333', zorder=5)
    assert not any(t.get_text() == 'B' for t in fig.texts)
    assert not any(t.get_text().startswith('Replication: KI') for t in fig.texts)
    assert not any(t.get_text().startswith('Bars: accepted estimates') for t in fig.texts)
    assert any(t.get_text() == 'Five channel transcripts' for t in fig.texts)
    assert any(t.get_text() == 'Four other transcripts' for t in fig.texts)
    assert not any(t.get_text() in HEADERS for t in fig.texts)
    assert ylabel.get_visible() and ylabel.get_text() == 'log2 fold change'
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = sorted({t.get_fontsize() for t in texts})
    assert min(sizes) >= 8
    fw, fh = fig.bbox.width, fig.bbox.height
    overflow = [t.get_text() for t in texts if (e := t.get_window_extent(renderer)).x0 < -.5
                or e.x1 > fw + .5 or e.y0 < -.5 or e.y1 > fh + .5]
    assert not overflow, overflow
    ext = [t.get_window_extent(renderer) for t in bottom_ax.get_xticklabels()]
    gaps_cm = [(ext[i + 1].x0 - ext[i].x1) / fig.dpi * 2.54 for i in range(len(ext) - 1)]
    assert min(gaps_cm) >= MIN_GAP_CM, ('x labels too close', gaps_cm)
    title_gaps_cm = {}
    for title, verdict, *_ in b.DATASETS:
        name = next(t for t in fig.texts if t.get_text() == DISPLAY.get(title, title))
        v = next(t for t in fig.texts if t.get_text() == row_labels[title] and t.get_position()[1] == name.get_position()[1])
        assert name.get_window_extent(renderer).x1 + 4 < v.get_window_extent(renderer).x0
        title_gaps_cm[DISPLAY.get(title, title)] = round((v.get_window_extent(renderer).x0 - name.get_window_extent(renderer).x1) / fig.dpi * 2.54, 3)
    # v07 checks: the asterisk count of every row, and each asterisk clears its bar, whisker, the zero line and the row
    # title by at least 2 px (300 dpi) and stays inside its axes.
    px300 = 300 / fig.dpi
    want = {t: len(marks.get(t, [])) for t in titles}
    got = {t: sum(1 for x in fig.axes[i].texts if x.get_text() == '*') for i, t in enumerate(titles)}
    assert got == want, (got, want)
    all_stars = [t for ax in fig.axes for t in ax.texts if t.get_text() == '*']
    zero_half = .8 / 2 / 72 * fig.dpi
    star_check = {}
    for (title, g), st in stars.items():
        ax = st.axes
        rows = table[table.dataset == title].set_index('gene')
        tt = next(x for x in fig.texts if x.get_text() == DISPLAY.get(title, title)).get_window_extent(renderer)
        j = b.GENES.index(g)
        e = st.get_window_extent(renderer)
        ab = ax.get_window_extent(renderer)
        cx, zero_y = ax.transData.transform((j, 0))
        bar_top = ax.transData.transform((j, max(0.0, float(rows.loc[g, 'log2FC']))))[1]
        whisker_top = ax.transData.transform((j, float(rows.loc[g, 'ci_high'])))[1]
        below = e.y0 - max(zero_y + zero_half, bar_top, whisker_top)
        above = tt.y0 - e.y1
        inside = bool(e.x0 >= ab.x0 and e.x1 <= ab.x1 and e.y0 >= ab.y0 and e.y1 <= ab.y1)
        star_check['%s | %s' % (title, g)] = dict(clear_below_px300=round(below * px300, 1), clear_title_px300=round(above * px300, 1),
                                                  centre_offset_px300=round(((e.x0 + e.x1) / 2 - cx) * px300, 2), inside_axes=inside)
        assert inside, (title, g, 'asterisk outside its axes')
        assert below * px300 >= MIN_CLEAR_PX300 and above * px300 >= MIN_CLEAR_PX300, (title, g, star_check['%s | %s' % (title, g)])
        assert abs((e.x0 + e.x1) / 2 - cx) * px300 < 1, (title, g)
    # v07: the note under the x labels, left-aligned with the panels: two lines; line 1 broken once if it does not fit.
    label_bottom = min(t.get_window_extent(renderer).y0 for t in bottom_ax.get_xticklabels())
    note_y = (label_bottom - NOTE_GAP_CM / 2.54 * fig.dpi) / fig.bbox.height
    note = fig.text(1.3 / b.WIDTH, note_y, NOTE, fontsize=8, ha='left', va='top', color='#444444')
    fig.canvas.draw()
    ne = note.get_window_extent(renderer)
    if ne.x1 > fig.bbox.width - .5:
        note.set_text(NOTE_TWO_LINES)
        fig.canvas.draw()
        ne = note.get_window_extent(renderer)
    assert ne.x0 >= -.5 and ne.x1 <= fig.bbox.width - .5 and ne.y0 >= .5 and ne.y1 <= fig.bbox.height, ('note does not fit', ne)
    note_gap_cm = (label_bottom - ne.y1) / fig.dpi * 2.54
    assert note_gap_cm >= MIN_NOTE_GAP_CM, ('note too close to the x labels', note_gap_cm)
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    assert min({t.get_fontsize() for t in texts}) >= 8
    overflow = [t.get_text() for t in texts if (e := t.get_window_extent(renderer)).x0 < -.5
                or e.x1 > fw + .5 or e.y0 < -.5 or e.y1 > fh + .5]
    assert not overflow, overflow
    OUT.mkdir(exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['svg.fonttype'] = 'none'
    for ext in ('png', 'svg', 'pdf'):
        meta = {'Software': None} if ext == 'png' else {'Date': None} if ext == 'svg' else {'CreationDate': None}
        fig.savefig(OUT / ('Figure_4.' + ext), dpi=300, metadata=meta)
    plt.close(fig)
    # Preserve the accepted table byte-for-byte, including its stored provenance labels.
    with (OUT / 'PLOTTED_BARS.tsv').open('xb') as f:
        f.write(accepted_path.read_bytes())
    assert b.sha(OUT / 'PLOTTED_BARS.tsv') == b.sha(accepted_path)
    check = dict(bars_checked=54, ci_pairs_checked=54, accepted_estimates_unchanged=True,
                 width_cm=b.WIDTH, height_cm=b.HEIGHT, minimum_font_pt=min(sizes), font_sizes_pt=sizes,
                 panel_letter=False, bulk_kcnc3_dagger_visible=True, edge_overflow=overflow,
                 top_line_title=False, explanatory_footnote=False,
                 colour_legends=['Five channel transcripts', 'Four other transcripts'],
                 x_tick_labels=[t.get_text() for t in bottom_ax.get_xticklabels()], x_tick_pt=TICK_PT,
                 min_x_label_gap_cm=round(min(gaps_cm), 3),
                 row_titles_drawn=[DISPLAY.get(t, t) for t, *_ in b.DATASETS], title_to_count_gap_cm=title_gaps_cm,
                 q_source=dict(path=str(Q_SOURCE), sha256=q_sha, column='q_core', population='Purkinje cell', age_weeks=12),
                 asterisks=star_check, n_asterisks=len(all_stars), asterisks_per_row=got, note_text=note.get_text(),
                 a1_limma_as_run_sha256=A1_LIMMA_AS_RUN, q_files=q_files, marked=marks,
                 note_lines=note.get_text().count('\n') + 1, note_gap_to_x_labels_cm=round(note_gap_cm, 3),
                 y_axis_label='log2 fold change',
                 source_sha256={k: dict(path=str(p), sha256=h) for k, (p, h) in b.SOURCES.items()},
                 accepted_bars_sha256=b.sha(accepted_path))
    with (OUT / 'Q_VALUES.tsv').open('x', encoding='utf-8') as f:
        cols = ['dataset', 'gene', 'log2FC', 'q', 'marked', 'n_genes_tested', 'source', 'source_sha256']
        f.write('\t'.join(cols) + '\n' + ''.join('\t'.join(str(r[c]) for c in cols) + '\n' for r in qrows))
    with (OUT / 'CONTENT_CHECK.json').open('x') as f:
        f.write(json.dumps(check, indent=2) + '\n')
    record = dict(task='Figure 4 v07b: asterisks for q < 0.05 in every single-nucleus row',
                  label_changes=['asterisks: ' + '; '.join('%s: %s' % (t, ', '.join(gs) or 'none') for t, gs in marks.items())] +
                                ['note under the x labels: ' + note.get_text().replace(chr(10), ' / '),
                                 'v05 labels and titles otherwise kept'],
                  script_sha256=b.sha(Path(__file__)), helper_sha256=SOURCE_SHA,
                  content_check=check, outputs={p.name: b.sha(p) for p in sorted(OUT.iterdir())},
                  finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    with (OUT / 'RECORD.json').open('x') as f:
        f.write(json.dumps(record, indent=2) + '\n')
    with (OUT / 'SHA256SUMS').open('x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS':
                f.write(f'{b.sha(p)}  {p.name}\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()

#%% debug
