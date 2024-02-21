from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from utils.io_utils import get_job_outputs


class Assembler:
    def __init__(self, cfg: DictConfig, vendor_dir: Path | str):
        self.cfg = cfg
        self.assembler_root = self._install(Path(vendor_dir))

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
    def run(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        self.pre_assembly_step(*args, **kwargs)
        self.run(**kwargs)
        self.post_assembly_step(*args, **kwargs)


def assembler_factory(assembler: str, cfg: DictConfig) -> Assembler:
    vendor_dir = ph.get_vendor_path()
    if assembler == "LJA":
        from asm.la_jolla import LaJolla

        return LaJolla(cfg=cfg, vendor_dir=vendor_dir)


@typechecked
def run(cfg: DictConfig, **kwargs) -> dict:
    exec_args = {
        "genome": None,
        "reads": None,
        "output_path": None,
    }

    exec_args.update(kwargs)

    assembler = assembler_factory(cfg.asm.name, cfg.asm)
    assembler(**exec_args)

    return get_job_outputs(exec_args)
