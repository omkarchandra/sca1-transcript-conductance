"""A1 step (a): download the three Purkinje-resolution datasets..

GSE215336 (Borgenheimer et al. 2022, Pcp2-ATXN1[82Q]): GSE215336_RAW.tar; GSE246183 human arm: GSE246183_230502_human_ctype.h5ad.gz;
GSE269430 (SCA7-266Q, Purkinje-enriched): GSE269430_RAW.tar and GSE269430_MULTI_seq_barcodes.csv.gz; plus every series matrix of
the three series. Each series into its own new folder real_v01/external/<GSE>/; 4 GB cap per file (Content-Length checked first);
URL, bytes and sha256 recorded. A tar is unpacked into its own new folder only if every member is a plain relative file. Nothing
downloaded is executed or unpickled here.

  python -I -B 02_replication/a1_purkinje/a1_download.py
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
import json
import re
import sys
import tarfile
import urllib.request

OUT = _ReleasePath(_release_paths.expand('${OUTPUT_ROOT}/a1_download'))
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/'
FILES = {'GSE215336': ['GSE215336_RAW.tar'],
         'GSE246183': ['GSE246183_230502_human_ctype.h5ad.gz'],
         'GSE269430': ['GSE269430_RAW.tar', 'GSE269430_MULTI_seq_barcodes.csv.gz']}
CAP = 4 * 1024 ** 3
UA = 'project/1.0 (python urllib)'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 24), b''):
            h.update(b)
    return h.hexdigest()


def req(url, method='GET'):
    return urllib.request.urlopen(urllib.request.Request(url, method=method, headers={'User-Agent': UA}), timeout=300)


def fetch(url, target):
    with req(url, 'HEAD') as r:
        length = int(r.headers.get('Content-Length', -1))
    if length > CAP:
        raise ValueError(f'{url} is {length} bytes, over the cap')
    h, n = hashlib.sha256(), 0
    with req(url) as r, open(target, 'xb') as f:
        for block in iter(lambda: r.read(1 << 22), b''):
            n += len(block)
            if n > CAP:
                raise ValueError(f'{url} grew past the cap')
            h.update(block)
            f.write(block)
    if length >= 0 and n != length:
        raise ValueError(f'{url}: {n} bytes, Content-Length {length}')
    return dict(url=url, file=target.name, bytes=n, sha256=h.hexdigest())


def unpack(path, dest):
    with tarfile.open(path) as tf:
        members = tf.getmembers()
        for m in members:
            parts = m.name.split('/')
            if not m.isfile() or m.name.startswith('/') or any(p in ('', '.', '..') for p in parts) or '\\' in m.name:
                raise ValueError(f'unsafe tar member {m.name!r} in {path.name}')
        dest.mkdir()
        for m in members:
            with tf.extractfile(m) as src, open(dest / m.name, 'xb') as dst:
                for block in iter(lambda: src.read(1 << 22), b''):
                    dst.write(block)
    return dict(archive=path.name, members=[m.name for m in members], folder=dest.name)


def main():
    if OUT.exists():
        raise ValueError('external/ already exists')
    OUT.parent.mkdir(exist_ok=True)
    OUT.mkdir()
    rec = dict(task='A1 step (a), downloads', script_sha256=sha(Path(__file__)), started_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
               python=sys.version.split()[0], series={})
    for gse, names in FILES.items():
        base = f'{GEO}{gse[:-3]}nnn/{gse}/'
        d = OUT / gse
        d.mkdir()
        with req(base + 'matrix/') as r:
            listing = r.read().decode('utf-8', 'replace')
        mats = sorted(set(re.findall(rf'href="({gse}[^"]*_series_matrix\.txt\.gz)"', listing)))
        rows = [fetch(base + 'matrix/' + m, d / m) for m in mats] + [fetch(base + 'suppl/' + n, d / n) for n in names]
        tars = [unpack(d / n, d / ('unpacked_' + n[:-4])) for n in names if n.endswith('.tar')]
        rec['series'][gse] = dict(base=base, files=rows, archives=tars)
        print(gse, [(x['file'], x['bytes']) for x in rows], [len(t['members']) for t in tars], flush=True)
    rec['finished_utc'] = datetime.now(timezone.utc).isoformat(timespec='seconds')
    with open(OUT / 'RECORD.json', 'x') as f:
        f.write(json.dumps(rec, indent=1) + '\n')
    files = sorted(p for p in OUT.rglob('*') if p.is_file())
    with open(OUT / 'SHA256SUMS', 'x') as f:
        for p in files:
            f.write(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n')
    print('SHA256SUMS', sha(OUT / 'SHA256SUMS'), len(files), 'files')


if __name__ == '__main__':
    main()

#%% debug
