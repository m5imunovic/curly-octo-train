from pathlib import Path
from typing import Dict, List, Optional, Union
from omegaconf import DictConfig
from typeguard import typechecked


@typechecked
def get_read_files(read_path: Path, suffix: Optional[List[str]] = None, override: bool = True) -> List[Path]:
    '''Get all files in the read_path with the given suffix.

    Args:
        read_path (Path):
            Path to the directory containing the reads. If the path is a file, it will be treated as a single read.
        suffix (Optional[list[str]], optional):
            List of suffixes to match. If None, files matching default suffixes will be returned.
        override (bool, optional):
            If True, only use the given suffixes. If False, use the default suffixes and append the given suffixes.

    Returns:
        list[Path]: sorted list of file paths in read_path with the given suffix.
    '''
    default_suffixes = ['.fasta', '.fa', '.fastq', '.fq']
    if not override:
        default_suffixes += suffix if suffix else []
    else:
        assert suffix is not None, 'Pattern must be specified when override is True'
        default_suffixes = suffix

    if read_path.is_file() and read_path.suffix in set(default_suffixes):
        return [read_path]
    elif read_path.is_dir():
        return sorted([p.resolve() for p in read_path.glob('**/*') if p.suffix in set(default_suffixes)])

    return []


@typechecked
def compose_cmd_params(params: Union[Dict, DictConfig]) -> str:
    '''Compose command line parameters from a dictionary.

    Args:
        params (dict):
            Params is a dictionary with following keys: `short`, `long`, `append`.
            'short' is a dictionary where key gets prepended with `-` and value gets appended to the command line.
            'long' is a dictionary where key gets prepended with `--` and value gets appended to the command line.
            'append' is a strings that gets appended to the command line.

    Returns:
        str: _description_
    '''
    short_params = ' '.join([f'-{k} {v}' for k, v in params['short'].items()]) if 'short' in params else ''
    long_params = ' '.join([f'--{k} {v}' for k, v in params['long'].items()]) if 'long' in params else ''
    append_params = f'{params["append"]}' if 'append' in params and params['append'] is not None else ''
    combined_params = ' '.join([short_params,  long_params,  append_params])
    return combined_params.strip()
