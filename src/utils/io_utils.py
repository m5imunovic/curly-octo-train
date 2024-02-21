import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Union

from omegaconf import DictConfig, ListConfig
from typeguard import typechecked


@typechecked
def get_read_files(
    read_path: Path, suffix: Optional[List[str]] = None, override: bool = True, regex: Optional[str] = None
) -> List[Path]:
    """Get all files in the read_path with the given suffix.

    Args:
        read_path (Path):
            Path to the directory containing the reads. If the path is a file, it will be treated as a single read.
        suffix (Optional[list[str]], optional):
            List of suffixes to match. If None, files matching default suffixes will be returned.
        override (bool, optional):
            If True, only use the given suffixes. If False, use the default suffixes and append the given suffixes.

    Returns:
        list[Path]: sorted list of file paths in read_path with the given suffix.
    """
    default_suffixes = [".fasta", ".fa", ".fastq", ".fq"]
    if not override:
        default_suffixes += suffix if suffix else []
    else:
        assert suffix is not None, "Pattern must be specified when override is True"
        default_suffixes = suffix

    results = []
    if read_path.is_file() and read_path.suffix in set(default_suffixes):
        results = [read_path]
    elif read_path.is_dir():
        results = sorted(p.absolute() for p in read_path.glob("**/*") if p.suffix in set(default_suffixes))

    if regex is not None:
        results = [p for p in results if re.search(regex, str(p))]

    return results


@typechecked
def compose_cmd_params(params: Union[Dict, DictConfig]) -> str:
    """Compose command line parameters from a dictionary.

    Args:
        params (dict):
            Params is a dictionary with following keys: `short`, `long`, `append`.
            'short' is a dictionary where key gets prepended with `-` and value gets appended to the command line.
            'long' is a dictionary where key gets prepended with `--` and value gets appended to the command line.
            'append' is a strings that gets appended to the command line.

    Returns:
        str: _description_
    """
    short_params = " ".join([f"-{k} {v}" for k, v in params["short"].items()]) if "short" in params else ""
    long_params = []
    if "long" in params:
        for k, v in params["long"].items():
            if isinstance(v, list) or isinstance(v, ListConfig):
                for item in v:
                    long_params.append(f"--{k} {item}")
            else:
                long_params.append(f"--{k} {v}")
    long_params = " ".join(long_params)
    append_params = f'{params["append"]}' if "append" in params and params["append"] is not None else ""
    combined_params = " ".join([short_params, long_params, append_params])
    return combined_params.strip()


class PathEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Path):
            return str(obj)
        return json.JSONEncoder.default(self, obj)


def get_job_outputs(exec_args) -> dict:
    assert "output_path" in exec_args, "output_path must be specified in exec_args"

    metadata_path = exec_args["output_path"] / "job_metadata"
    metadata_path.mkdir(parents=True, exist_ok=True)

    with open(metadata_path / "exec_args.json", "w") as f:
        json.dump(exec_args, f, indent=2, cls=PathEncoder)

    produced_files = {
        "metadata": list(metadata_path.glob("**/*")),
        "artifacts": list(f for f in exec_args["output_path"].glob("**/*") if f.is_file()),
    }

    return produced_files
