from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig
from typeguard import typechecked


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
    def run(self, reads_path: Path, save_path: Path, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        # TODO: refactor reads simulator to use the same pattern
        self.pre_assembly_step(*args, **kwargs)
        self.run(**kwargs)
        self.post_assembly_step(*args, **kwargs)


def assembler_factory(assembler: str, cfg: DictConfig) -> Assembler:
    vendor_dir = cfg.paths.vendor_dir
    if assembler == "LJA":
        from asm.la_jolla import LaJolla

        return LaJolla(cfg=cfg.asm, vendor_dir=vendor_dir)


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        "genome": None,
        "reads": None,
        "out_path": None,
    }

    exec_args.update(kwargs)

    assembler = assembler_factory(cfg.asm.name, cfg)
    assembler(**exec_args)
