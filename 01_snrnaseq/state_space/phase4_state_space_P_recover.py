"""Report-only recovery of the single failed v05g P output bundle.

This script copies already-produced report/matrix artifacts. It never reruns
scoring and never writes to the failed source directory. A real invocation
requires Python -I -B, a review-byte-bound GO file, and <AUTHOR>'s exact
"yes recover P" relay. Synthetic fixtures call helpers with synthetic paths;
they never call main().

GO schema:
  scope: P_report_recovery
  script_sha256: {phase4_state_space_P_recover.py: <full sha256>}
  byte_clearance: yes
  author_verbatim: "yes recover P"
  relayed_by and relayed_at_cdt: non-empty provenance strings
"""
from __future__ import annotations

# Load release-local path configuration without relying on the working directory.
import importlib.util as _release_importlib
from pathlib import Path as _ReleasePath
_release_config = next(p / "config" / "paths.py" for p in _ReleasePath(__file__).resolve().parents if (p / "config" / "paths.py").is_file())
_release_spec = _release_importlib.spec_from_file_location("release_paths", _release_config)
_release_paths = _release_importlib.module_from_spec(_release_spec)
_release_spec.loader.exec_module(_release_paths)


import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional


WORKSPACE = Path(_release_paths.expand('${WORK_ROOT}/state_space'))
SOURCE_DIR = WORKSPACE / 'interface/P_v05g_01'
RECOVERY_DIR = WORKSPACE / 'interface/P_v05g_01_recovery01'
CONDA = Path(_release_paths.expand('${CONDA_BIN}'))
ENV = WORKSPACE / 'env/single_cell_ucell'
JOB_ID = 9286379
SCRIPT_PATH = Path(__file__).resolve()

POPULATION_SLUGS = ('Purkinje_cell', 'Bergmann_glia')
TRIPLET_FILES = tuple(
    name
    for slug in POPULATION_SLUGS
    for name in (
        f'state_P_{slug}.npy',
        f'state_P_{slug}_nuclei.tsv',
        f'state_P_{slug}_features.tsv',
    )
)
COPY_FILES = (*TRIPLET_FILES, 'RECORD.json')
FAILED_BUILDER_NAME = 'phase4_state_space_P_v05g.py'
FAILED_BUILDER_SHA256 = '6d995c95402abd33a5c1b08711e99b4c92ea3886b7b0ac740003428d5222d4e2'
FAILED_RECORD_SHA256 = '123140b32717e97d6d63d63d6a479a79293585695ebd99d133b8f759ac8392da'
FAILURE_CAUSE = 'TypeError: dict() got multiple values for keyword argument \'T6b_script_sha256\''
TRIPLET_SHA256 = {
    'state_P_Purkinje_cell.npy': 'e5dd62b8ab8305e23211408a151fd6cf39ee78bf1dc14bd25d8aec4a64e6c2ea',
    'state_P_Purkinje_cell_nuclei.tsv': '9882ddc0dff85c55f1daea499c1791801ad87a3864b25df9c71ad64ebabaf82a',
    'state_P_Purkinje_cell_features.tsv': 'd141a67e7b4f462da1136164d43aea8cd5a20c0c18f4a472d6e920daea8047b7',
    'state_P_Bergmann_glia.npy': 'd54c151635f2ea9312687fd473103b96c711f776948996b34510c30f04dd9de3',
    'state_P_Bergmann_glia_nuclei.tsv': '3cb5d01f46fafe8341bc655828734ea1d0fd25a385b31247b5215f27447b2e2a',
    'state_P_Bergmann_glia_features.tsv': 'd141a67e7b4f462da1136164d43aea8cd5a20c0c18f4a472d6e920daea8047b7',
}

# Explicitly named so synthetic fixtures can neutralize the real bindings
# immediately after importing this module.
CFG = {
    'source_dir': SOURCE_DIR,
    'recovery_dir': RECOVERY_DIR,
    'conda': CONDA,
    'env': ENV,
    'failed_job_id': JOB_ID,
}
PINS = {
    'failed_builder_name': FAILED_BUILDER_NAME,
    'failed_builder_sha256': FAILED_BUILDER_SHA256,
    'record_sha256': FAILED_RECORD_SHA256,
    'triplet_sha256': TRIPLET_SHA256,
}


class RecoveryRefusal(RuntimeError):
    """A read-only preflight or GO gate did not match the declared pins."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(16 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json_exclusive(path: Path, value: dict[str, Any]) -> None:
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n'
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(payload)


def write_bytes_exclusive(path: Path, payload: bytes) -> None:
    with Path(path).open('xb') as stream:
        stream.write(payload)


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError) as exc:
        raise RecoveryRefusal(f'{label} is unreadable or invalid JSON: {exc}') from exc
    if not isinstance(value, dict):
        raise RecoveryRefusal(f'{label} must contain a JSON object')
    return value


def check_go(go_file: Path, script_path: Path = SCRIPT_PATH) -> dict[str, Any]:
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise RecoveryRefusal('refusing: assertion checks require non-optimized Python')
    if sys.flags.isolated != 1 or not sys.dont_write_bytecode:
        raise RecoveryRefusal('refusing unless launched with python -I -B')
    gate = _load_json(go_file, 'GO file')
    expected = {Path(script_path).name: sha256_file(Path(script_path))}
    if gate.get('scope') != 'P_report_recovery':
        raise RecoveryRefusal('GO scope must be P_report_recovery')
    if gate.get('script_sha256') != expected:
        raise RecoveryRefusal(f'GO must bind this exact recovery script hash: {expected}')
    if gate.get('byte_clearance') != 'yes':
        raise RecoveryRefusal('byte review is absent')
    if gate.get('author_verbatim') != 'yes recover P':
        raise RecoveryRefusal('<AUTHOR> verbatim approval must be exactly "yes recover P"')
    if not gate.get('relayed_by') or not gate.get('relayed_at_cdt'):
        raise RecoveryRefusal('GO relay provenance is absent')
    return gate


def _require_regular_file(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file():
        raise RecoveryRefusal(f'{label} must be a regular non-symlink file: {path}')


def validate_source(source_dir: Path, recovery_dir: Path,
                    pins: dict[str, Any]) -> dict[str, Any]:
    source = Path(source_dir)
    recovery = Path(recovery_dir)
    if source.is_symlink() or not source.is_dir():
        raise RecoveryRefusal(f'failed source must be a real directory: {source}')
    if recovery.exists() or recovery.is_symlink():
        raise RecoveryRefusal(f'recovery directory already exists: {recovery}')
    stop_path = source / 'STOP.json'
    if stop_path.exists() or stop_path.is_symlink():
        raise RecoveryRefusal(f'STOP.json is present: {stop_path}')

    failed_path = source / 'FAILED.json'
    record_path = source / 'RECORD.json'
    _require_regular_file(failed_path, 'FAILED.json')
    _require_regular_file(record_path, 'RECORD.json')
    failed = _load_json(failed_path, 'FAILED.json')
    traceback_text = failed.get('traceback')
    cause = FAILURE_CAUSE
    if not isinstance(traceback_text, str) or cause not in traceback_text:
        raise RecoveryRefusal('FAILED.json does not contain the pinned duplicate-key TypeError')
    failed_scripts = failed.get('script_sha256')
    expected_builder = pins['failed_builder_sha256']
    builder_name = pins['failed_builder_name']
    if not isinstance(failed_scripts, dict) or failed_scripts.get(builder_name) != expected_builder:
        raise RecoveryRefusal('FAILED.json does not bind the approved v05g builder hash')

    record_hash = sha256_file(record_path)
    if record_hash != pins['record_sha256']:
        raise RecoveryRefusal(f'RECORD.json hash mismatch: {record_hash}')
    record = _load_json(record_path, 'RECORD.json')
    if record.get('mode') != 'real' or record.get('scope') != 'gene_set_space_P':
        raise RecoveryRefusal('RECORD.json mode/scope do not match the failed P run')
    if record.get('completed') is not True or record.get('analysis_disposition') != 'PASS':
        raise RecoveryRefusal('RECORD.json is not a completed PASS analysis record')
    if not isinstance(record.get('QC6_Bergmann'), dict):
        raise RecoveryRefusal('RECORD.json lacks the QC6_Bergmann object to preserve')

    triplet_hashes: dict[str, str] = {}
    expected_triplets = pins['triplet_sha256']
    if set(expected_triplets) != set(TRIPLET_FILES):
        raise RecoveryRefusal('internal triplet pins do not name exactly the six interface files')
    for name in TRIPLET_FILES:
        source_file = source / name
        _require_regular_file(source_file, name)
        actual = sha256_file(source_file)
        if actual != expected_triplets[name]:
            raise RecoveryRefusal(f'{name} hash mismatch: {actual}')
        triplet_hashes[name] = actual

    return {
        'failed': failed,
        'failed_sha256': sha256_file(failed_path),
        'record': record,
        'record_sha256': record_hash,
        'triplet_sha256': triplet_hashes,
    }


def _source_snapshot(source_dir: Path) -> tuple[tuple[tuple[Any, ...], ...], dict[str, str]]:
    source = Path(source_dir)
    entries = []
    for child in sorted(source.iterdir(), key=lambda item: item.name):
        st = child.lstat()
        entries.append((child.name, st.st_mode, st.st_size, st.st_mtime_ns))
    names = ('FAILED.json', *COPY_FILES)
    hashes = {name: sha256_file(source / name) for name in names}
    return tuple(entries), hashes


def _copy_exclusive(source: Path, destination: Path, expected_sha256: str) -> None:
    digest = hashlib.sha256()
    with Path(source).open('rb') as src, Path(destination).open('xb') as dst:
        for block in iter(lambda: src.read(16 << 20), b''):
            digest.update(block)
            dst.write(block)
    copied_hash = digest.hexdigest()
    if copied_hash != expected_sha256 or sha256_file(destination) != expected_sha256:
        raise RecoveryRefusal(f'copy hash verification failed for {destination.name}')


def perform_recovery(source_dir: Path, recovery_dir: Path, cfg: dict[str, Any],
                     pins: dict[str, Any], script_path: Path,
                     explicit_bytes: Optional[bytes] = None) -> dict[str, Any]:
    """Copy the pinned failed-run bundle and write provenance; never score data."""
    source = Path(source_dir)
    recovery = Path(recovery_dir)
    validated = validate_source(source, recovery, pins)
    before = _source_snapshot(source)

    if explicit_bytes is None:
        command = [str(cfg['conda']), 'list', '--explicit', '-p', str(cfg['env'])]
        completed = subprocess.run(command, check=True, capture_output=True)
        explicit_bytes = completed.stdout
    if not explicit_bytes:
        raise RecoveryRefusal('builder-equivalent conda explicit output is empty')

    # Repeat all pins and the source snapshot immediately before creating output.
    validated_again = validate_source(source, recovery, pins)
    if validated_again['failed_sha256'] != validated['failed_sha256']:
        raise RecoveryRefusal('FAILED.json changed during preflight')
    if validated_again['record_sha256'] != validated['record_sha256']:
        raise RecoveryRefusal('RECORD.json changed during preflight')
    if validated_again['triplet_sha256'] != validated['triplet_sha256']:
        raise RecoveryRefusal('triplet files changed during preflight')
    if _source_snapshot(source) != before:
        raise RecoveryRefusal('failed source directory changed during preflight')

    recovery.mkdir(exist_ok=False)
    copied_hashes: dict[str, str] = {}
    for name in COPY_FILES:
        expected = (validated['record_sha256'] if name == 'RECORD.json'
                    else validated['triplet_sha256'][name])
        _copy_exclusive(source / name, recovery / name, expected)
        copied_hashes[name] = expected

    if _source_snapshot(source) != before:
        raise RecoveryRefusal('failed source directory changed while copying; preserve partial recovery')
    for name, expected in copied_hashes.items():
        if sha256_file(recovery / name) != expected:
            raise RecoveryRefusal(f'recovery copy changed after verification: {name}')

    result = {
        'mode': validated['record']['mode'],
        'scope': validated['record']['scope'],
        'completed': True,
        'analysis_disposition': validated['record']['analysis_disposition'],
        'recovered_from_failed_run': True,
        'failed_job_id': int(cfg['failed_job_id']),
        'source_directory': str(source),
        'failed_json_sha256': validated['failed_sha256'],
        'record_sha256': validated['record_sha256'],
        'triplet_sha256': validated['triplet_sha256'],
        'recovery_script_sha256': sha256_file(script_path),
        'cause': FAILURE_CAUSE,
        'QC6_Bergmann': validated['record']['QC6_Bergmann'],
    }
    write_json_exclusive(recovery / 'RESULT.json', result)
    write_bytes_exclusive(recovery / 'CONDA_EXPLICIT.txt', explicit_bytes)
    manifest = ''.join(
        f"{validated['triplet_sha256'][name]}  {name}\n" for name in TRIPLET_FILES
    )
    with (recovery / 'SHA256SUMS.txt').open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(manifest)
    return result


def main(argv: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--go-file', type=Path, required=True)
    args = parser.parse_args(argv)
    check_go(args.go_file)
    result = perform_recovery(
        CFG['source_dir'], CFG['recovery_dir'], CFG, PINS, SCRIPT_PATH
    )
    print(json.dumps({
        'status': 'PASS',
        'recovery_directory': str(CFG['recovery_dir']),
        'recovery_script_sha256': result['recovery_script_sha256'],
        'triplet_sha256': result['triplet_sha256'],
    }, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
