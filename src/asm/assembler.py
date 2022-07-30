from abc import abstractmethod
from pathlib import Path
from typing import Optional

import hydra
# from omegaconf import OmegaConf
from typeguard import typechecked

from utils import path_helpers as ph


# OmegaConf.register_new_resolver('project_root', ph.project_root_append)

class Assembler:
    def __init__(self, cfg: dict, vendor_dir: Path):
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


def assembler_factory(assembler: str, params: Optional[dict] = None) -> Assembler:
    vendor_dir: Path = ph.get_vendor_path()
    if assembler == 'LJA':
        from asm.la_jolla import LaJolla
        return LaJolla(cfg=params, vendor_dir=vendor_dir)


def run(cfg):
    assembler_cfg = dict(cfg)
    assembler = assembler_factory(assembler_cfg['name'], assembler_cfg)
    kwargs = {
        'reads_path': ph.get_simulated_data_path() / assembler_cfg['species'],
        'out_path': ph.get_assemblies_path() / assembler_cfg['experiment']
    }

    assembler(**kwargs)


@hydra.main(version_base=None, config_path='../../config/asm', config_name='la_jolla')
def main(cfg):
    print("Running assembler step...")
    run(cfg)


if __name__ == "__main__":
    main()