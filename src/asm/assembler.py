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


def assembler_factory(assembler: str, cfg: DictConfig) -> Assembler:
    vendor_dir: Path = cfg.paths.vendor_dir
    if assembler == 'LJA':
        from asm.la_jolla import LaJolla
        return LaJolla(cfg=cfg.asm, vendor_dir=vendor_dir)


def assembly_experiment_path(cfg: DictConfig) -> Path:
    return cfg.paths.assemblies_dir / cfg.asm.experiment


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        'reads_path': Path(cfg.paths.simulated_data_dir) / cfg.species_name,
        'out_path': assembly_experiment_path(cfg)
    }

    exec_args.update(kwargs)

    assembler = assembler_factory(cfg.asm.name, cfg)
    assembler(**exec_args)


@hydra.main(version_base=None, config_path='../../config/asm', config_name='la_jolla')
def main(cfg):
    print("Running assembler step...")

    run(cfg)


if __name__ == "__main__":
    main()