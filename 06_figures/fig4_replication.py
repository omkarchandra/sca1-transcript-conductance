"""Plot accepted nine-gene replication estimates and CIs as standalone Figure 4."""

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)

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
SOURCE = HERE.parent / _release_paths.expand('${REPO_ROOT}/06_figures/replication_values.py')
SOURCE_SHA = '68d63352f564fba72c53431733eae7c1d2d8e74c4e262639fa97b77f850f33a5'
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/fig4_replication'))


def main():
    assert not OUT.exists(), 'exclusive-create output; do not overwrite'
    spec = importlib.util.spec_from_file_location('accepted_replication', SOURCE)
    b = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(b)
    assert b.sha(SOURCE) == SOURCE_SHA
    table = b.load()
    accepted_path = SOURCE.parent / _release_paths.expand('${REPO_ROOT}/results/replication/PLOTTED_BARS.tsv')
    accepted = b.pd.read_csv(accepted_path, sep='\t')
    assert list(table.gene) == list(accepted.gene) and list(table.dataset) == list(accepted.dataset)
    for col in ('log2FC', 'ci_low', 'ci_high'):
        assert np.allclose(table[col], accepted[col], rtol=0, atol=1e-14)
    fig = plt.figure(figsize=(b.WIDTH / 2.54, b.HEIGHT / 2.54))
    b.draw(fig, table)
    # Only the panel letter is removed; all data geometry stays as accepted.
    for text in list(fig.texts):
        if text.get_text() == 'B':
            text.remove()
    bulk = fig.axes[4]
    marker = bulk.text(1, .35, '†', fontsize=9, ha='center', va='bottom', color='#333333', zorder=5)
    assert marker.get_visible() and marker.get_text() == '†'
    assert not any(t.get_text() == 'B' for t in fig.texts)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    sizes = sorted({t.get_fontsize() for t in texts})
    assert min(sizes) >= 8
    fw, fh = fig.bbox.width, fig.bbox.height
    overflow = [t.get_text() for t in texts if (e := t.get_window_extent(renderer)).x0 < -.5
                or e.x1 > fw + .5 or e.y0 < -.5 or e.y1 > fh + .5]
    assert not overflow, overflow
    for title, verdict, *_ in b.DATASETS:
        name = next(t for t in fig.texts if t.get_text() == title)
        v = next(t for t in fig.texts if t.get_text() == verdict and t.get_position()[1] == name.get_position()[1])
        assert name.get_window_extent(renderer).x1 + 4 < v.get_window_extent(renderer).x0
    OUT.mkdir(exist_ok=False)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['svg.fonttype'] = 'none'
    for ext in ('png', 'svg', 'pdf'):
        meta = {'Software': None} if ext == 'png' else {'Date': None} if ext == 'svg' else {'CreationDate': None}
        fig.savefig(OUT / ('Figure_4.' + ext), dpi=300, metadata=meta)
    plt.close(fig)
    table.to_csv(OUT / 'PLOTTED_BARS.tsv', sep='\t', index=False)
    check = dict(bars_checked=54, ci_pairs_checked=54, accepted_estimates_unchanged=True,
                 width_cm=b.WIDTH, height_cm=b.HEIGHT, minimum_font_pt=min(sizes), font_sizes_pt=sizes,
                 panel_letter=False, bulk_kcnc3_dagger_visible=True, edge_overflow=overflow,
                 source_sha256={k: dict(path=str(p), sha256=h) for k, (p, h) in b.SOURCES.items()},
                 accepted_bars_sha256=b.sha(accepted_path))
    (OUT / 'CONTENT_CHECK.json').write_text(json.dumps(check, indent=2) + '\n')
    record = dict(script_sha256=b.sha(Path(__file__)), helper_sha256=SOURCE_SHA,
                  content_check=check, outputs={p.name: b.sha(p) for p in sorted(OUT.iterdir())},
                  finished_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'))
    (OUT / 'RECORD.json').write_text(json.dumps(record, indent=2) + '\n')
    with (OUT / 'SHA256SUMS').open('x') as f:
        for p in sorted(OUT.iterdir()):
            if p.name != 'SHA256SUMS':
                f.write(f'{b.sha(p)}  {p.name}\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()

#%% debug
