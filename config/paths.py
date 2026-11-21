"""Resolve repository paths and explicit local YAML configuration placeholders."""
from functools import lru_cache
from pathlib import Path
import argparse
import os
import re
import shlex
import yaml

KEYS = ('WORK_ROOT',
        'DATA_ROOT', 'OUTPUT_ROOT', 'CONDA_BIN', 'CONDA_ENV',
        'MODEL_ROOT', 'PYTHON_NEURON')


def repo_root():
    return Path(__file__).resolve().parent.parent


@lru_cache(maxsize=None)
def load(config_file=None):
    file = Path(config_file or os.environ.get('LOCAL_PATH_CONFIG')
                or Path(__file__).with_name('local.yaml')).expanduser().resolve()
    with file.open(encoding='utf-8') as f:
        values = yaml.safe_load(f)
    if not isinstance(values, dict):
        raise ValueError('Local path configuration must be a mapping')
    missing = [key for key in KEYS if not isinstance(values.get(key), str) or not values[key].strip()]
    if missing:
        raise ValueError('Missing local configuration: ' + ', '.join(missing))
    return values


def value(key, config_file=None):
    if key == 'REPO_ROOT':
        return str(repo_root())
    return os.environ.get(key) or load(config_file)[key]


def path(key, *parts, config_file=None):
    return Path(value(key, config_file)).expanduser().joinpath(*parts)


def expand(text, config_file=None):
    return re.sub(r'\$\{([A-Z_]+)\}', lambda match: value(match[1], config_file), text)


def shell_exports(config_file=None):
    file = Path(config_file or os.environ.get('LOCAL_PATH_CONFIG')
                or Path(__file__).with_name('local.yaml')).expanduser().resolve()
    values = load(str(file))
    values = {key: os.environ.get(key) or values[key] for key in KEYS}
    values.update(LOCAL_PATH_CONFIG=str(file), RELEASE_ROOT=str(repo_root()), REPO_ROOT=str(repo_root()))
    return '\n'.join('export ' + key + '=' + shlex.quote(str(value)) for key, value in values.items())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shell', action='store_true', required=True)
    parser.add_argument('--config')
    args = parser.parse_args()
    print(shell_exports(args.config))

