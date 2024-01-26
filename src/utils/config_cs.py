from dataclasses import dataclass
from typing import Any

from hydra.core.config_store import ConfigStore

from utils import path_cs


@dataclass
class Config:
    hydra: Any
    paths: path_cs.PathsConfig
    species_name: dict
    seed: int
    tags: str | None
    experiment: str | None
    date_mm_dd: str | None

    reference: dict
    reads: dict
    asm: dict
    graph: dict


def register_configs():
    cs = ConfigStore.instance()
    cs.store(name="base_reference", node=Config)
    path_cs.register_config()
