from abc import abstractmethod
from pathlib import Path
from typing import Optional

import hydra
from omegaconf import DictConfig
from typeguard import typechecked

from utils import path_helpers as ph


class Assembler:
    def __init__(self, cfg: DictConfig, vendor_dir: Path):
        self.cfg = cfg
        self.assembler_root = self._install(vendor_dir)

    @abstractmethod
    @typechecked
    def _install(self, vendor_dir: Path):
        pass

    def pre_assembly_step(self, *args, **kwargs):
        pass

    def post_assembly_step(self, *args, **kwargs):
        pass

    @abstractmethod
    @typechecked
    def run(self, reads_path: Path, save_path: Path, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        self.pre_assembly_step(*args, **kwargs)
        self.run(**kwargs)
        self.post_assembly_step(*args, **kwargs)


def assembler_factory(assembler: str, params: DictConfig) -> Assembler:
    vendor_dir: Path = ph.get_vendor_path()
    if assembler == 'LJA':
        from asm.la_jolla import LaJolla
        return LaJolla(cfg=params, vendor_dir=vendor_dir)


def assembly_experiment_path(cfg: DictConfig) -> Path:
    return ph.get_assemblies_path() / cfg.experiment


def run(cfg: DictConfig, **kwargs):
    assembler = assembler_factory(cfg['name'], cfg)
    assembler(**kwargs)


@hydra.main(version_base=None, config_path='../../config/asm', config_name='la_jolla')
def main(cfg):
    print("Running assembler step...")

    kwargs = {
        'reads_path': ph.get_simulated_data_path() / cfg.species,
        'out_path': assembly_experiment_path()
    }

    run(cfg, **kwargs)


if __name__ == "__main__":
    main()