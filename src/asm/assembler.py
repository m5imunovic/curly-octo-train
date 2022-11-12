from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig, OmegaConf
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
    return Path(cfg.paths.assemblies_dir) / cfg.experiment


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        'reads_path': Path(cfg.paths.reads_dir) / cfg.species_name,
        'ref_path': Path(cfg.paths.ref_dir) / cfg.species_name / 'chromosomes',
        'out_path': assembly_experiment_path(cfg),
        'suffix': OmegaConf.to_container(cfg.suffix)
    }

    exec_args.update(kwargs)

    assembler = assembler_factory(cfg.asm.name, cfg)
    assembler(**exec_args)
