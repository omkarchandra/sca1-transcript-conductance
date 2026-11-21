"""P2 step (b): GEO downloads..

v02: v01 (9cedd53d) stopped in series_design on a blank line of the GSE75778 series matrix (IndexError); v02 skips blank
lines and writes to a new real_v01/external_v02/. The partial v01 folder real_v01/external/ is kept as it is.

Refuses to start unless the control set (step a) is on disk at its pinned hash. For GSE75778 and GSE122099: the GEO listing of
matrix/ and suppl/, every series matrix and every supplementary file (2 GB cap per file, checked from Content-Length before
the body is fetched), each written once into its own new folder real_v01/external_v02/<GSE>/files/ with URL, bytes and sha256.
A tar archive is listed first and unpacked into its own new folder only if every member is a plain relative file. From each
series matrix: the series title, summary and design lines and one row per sample with its title, source, characteristics and
supplementary files (DESIGN.tsv). Downloaded files are only read as text or listed; nothing in them is executed.

  python -I -B 02_replication/p2_bulk/p2_download.py
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
from html.parser import HTMLParser
import csv
import gzip
import hashlib
import json
import re
import sys
import tarfile
import urllib.error
import urllib.parse
import urllib.request

HERE = Path(__file__).resolve().parent
OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/p2_download'))
CONTROL = HERE / 'real_v01' / 'control'
CONTROL_PINS = {'CONTROL_SET.tsv': '1968a1fdfbcb51371eb6d3cf1f725047d36fceeb790abd78fc412c5504bcc289',
                'SHA256SUMS': '5f2043fee71304323223b47a9d76e4771ee3847a21ee6e8768fb79e19bcec021'}
SERIES = ('GSE75778', 'GSE122099')
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
CAP = 2 * 1024 ** 3
USER_AGENT = 'project/1.0 (python urllib)'
NAME_OK = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._+-]*$')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.hrefs += [v for k, v in attrs if k == 'href' and v]


def request(url, method='GET'):
    return urllib.request.urlopen(urllib.request.Request(url, method=method, headers={'User-Agent': USER_AGENT}), timeout=120)


def listing(url, save_as):
    try:
        with request(url) as r:
            body = r.read(5 * 1024 ** 2)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None                     # no such GEO folder; recorded, not an error
        raise
    with open(save_as, 'xb') as f:
        f.write(body)
    p = Links()
    p.feed(body.decode('utf-8', 'replace'))
    names = []
    for h in p.hrefs:
        if h.startswith(('?', '/', '#', 'http', '..')) or h.endswith('/'):
            continue
        name = urllib.parse.unquote(h)
        if not NAME_OK.match(name):
            raise ValueError(f'unexpected file name in {url}: {name!r}')
        names.append(name)
    return sorted(set(names))


def fetch(url, target):
    with request(url, 'HEAD') as r:
        length = int(r.headers.get('Content-Length', -1))
        head = dict(r.headers.items())
    if length > CAP:
        return dict(url=url, file=target.name, status='skipped_over_cap', bytes=length, sha256='', last_modified=head.get('Last-Modified', ''))
    h, n = hashlib.sha256(), 0
    with request(url) as r, open(target, 'xb') as f:
        for block in iter(lambda: r.read(1 << 20), b''):
            n += len(block)
            if n > CAP:
                raise ValueError(f'{url} grew past the cap')
            h.update(block)
            f.write(block)
    if length >= 0 and n != length:
        raise ValueError(f'{url}: {n} bytes, Content-Length {length}')
    return dict(url=url, file=target.name, status='downloaded', bytes=n, sha256=h.hexdigest(), last_modified=head.get('Last-Modified', ''))


def unpack_tar(path, dest):
    with tarfile.open(path) as tf:
        members = tf.getmembers()
        for m in members:
            parts = m.name.split('/')
            if not (m.isfile() or m.isdir()) or m.name.startswith('/') or any(p in ('', '.', '..') for p in parts) or '\\' in m.name:
                return dict(archive=path.name, unpacked=False, reason=f'unsafe member {m.name!r}', members=len(members))
        dest.mkdir()
        for m in members:
            target = dest / m.name
            if m.isdir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with tf.extractfile(m) as src, open(target, 'xb') as dst:
                dst.write(src.read())
    return dict(archive=path.name, unpacked=True, members=len(members), folder=dest.name)


def series_design(path):
    series, samples = {}, {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                break
            if not line.strip():
                continue
            row = next(csv.reader([line.rstrip('\n')], delimiter='\t'))
            key, vals = row[0], row[1:]
            if key.startswith('!Series_') and key in ('!Series_title', '!Series_summary', '!Series_overall_design', '!Series_supplementary_file'):
                series.setdefault(key[8:], []).extend(vals)
            elif key.startswith('!Sample_'):
                field = key[8:]
                k = 1
                while f'{field}_{k}' in samples:
                    k += 1
                samples[f'{field}_{k}' if field in samples or field.startswith('characteristics') else field] = vals
    return series, samples


def main():
    for name, want in CONTROL_PINS.items():
        if sha(CONTROL / name) != want:
            raise ValueError(f'control set not at its pinned hash: {name}')
    if OUT.exists():
        raise ValueError('external_v02/ already exists')
    started = datetime.now(timezone.utc)
    OUT.mkdir(parents=False, exist_ok=False)
    record = dict(task='P2 step (b), GEO downloads', script_sha256=sha(Path(__file__)), control_pins=CONTROL_PINS,
                  started_utc=started.isoformat(timespec='seconds'), python=sys.version.split()[0], user_agent=USER_AGENT, series={})
    for gse in SERIES:
        base = f'{GEO}{gse[:-3]}nnn/{gse}/'
        d = OUT / gse
        d.mkdir()
        (d / 'files').mkdir()
        rows, design = [], {}
        for sub in ('matrix', 'suppl'):
            names = listing(base + sub + '/', d / f'listing_{sub}.html')
            if names is None:
                rows.append(dict(folder=sub, url=base + sub + '/', file='', status='folder_not_found', bytes=0, sha256='', last_modified=''))
                continue
            for name in names:
                rows.append(dict(folder=sub, **fetch(base + sub + '/' + urllib.parse.quote(name), d / 'files' / name)))
        archives = []
        for r in rows:
            if r['status'] == 'downloaded' and re.search(r'\.tar(\.gz)?$', r['file']):
                archives.append(unpack_tar(d / 'files' / r['file'], d / ('unpacked_' + r['file'].split('.tar')[0])))
        for r in rows:
            if r['folder'] == 'matrix' and r['status'] == 'downloaded' and r['file'].endswith('_series_matrix.txt.gz'):
                series, samples = series_design(d / 'files' / r['file'])
                design[r['file']] = dict(series=series, n_samples=len(samples.get('geo_accession', [])))
                n = len(samples.get('geo_accession', []))
                with open(d / f'DESIGN_{r["file"].replace("_series_matrix.txt.gz", "")}.tsv', 'x', newline='') as f:
                    w = csv.writer(f, delimiter='\t', lineterminator='\n')
                    keys = [k for k in samples if len(samples[k]) == n]
                    w.writerow(keys)
                    for i in range(n):
                        w.writerow([samples[k][i] for k in keys])
        with open(d / 'DOWNLOADS.tsv', 'x', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['folder', 'url', 'file', 'status', 'bytes', 'sha256', 'last_modified'], delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
        record['series'][gse] = dict(base=base, files=rows, archives=archives, design=design)
    record['finished_utc'] = datetime.now(timezone.utc).isoformat(timespec='seconds')
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(record, indent=1) + '\n')
    files = sorted(p for p in OUT.rglob('*') if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n')
    for gse, s in record['series'].items():
        print(gse, [(r['file'], r['status'], r['bytes']) for r in s['files']], s['archives'])
        for name, des in s['design'].items():
            print(' ', name, des['n_samples'], {k: v[:2] for k, v in des['series'].items() if k != 'summary'})
    print('SHA256SUMS', sha(OUT / 'SHA256SUMS'), len(files), 'files')


if __name__ == '__main__':
    main()

#%% debug
