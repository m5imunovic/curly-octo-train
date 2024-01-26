from dataclasses import dataclass
from pathlib import Path

from hydra.core.config_store import ConfigStore


@dataclass
class PathsConfig:
    root_dir: Path
    vendor_dir: Path
    data_dir: Path
    ref_dir: Path
    reads_dir: Path
    assemblies_dir: Path
    datasets_dir: Path
    output_dir: Path
    log_dir: Path


def register_config():
    cs = ConfigStore.instance()
    cs.store(name="base_paths", group="paths", node=PathsConfig)
