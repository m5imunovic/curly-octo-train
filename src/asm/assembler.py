import logging
from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from utils.io_utils import get_job_outputs

logger = logging.getLogger(__name__)
MAX_RETRIES = 2


class Assembler:
    def __init__(self, cfg: DictConfig, vendor_dir: Path | str):
        self.cfg = cfg
        self.assembler_root = self._install(Path(vendor_dir))
        self.max_retries = MAX_RETRIES if "max_retries" not in self.cfg else self.cfg.max_retries
        self.read_suff = [".fastq", ".fq"] if "read_suffix_filter" not in self.cfg else self.cfg.read_suffix_filter

    @abstractmethod
    @typechecked
    def _install(self, vendor_dir: Path):
        pass

    def pre_assembly_step(self, *args, **kwargs):
        pass

    def post_assembly_step(self, *args, **kwargs):
        pass

    @property
    def max_retries(self):
        return self._max_retries

    @max_retries.setter
    def max_retries(self, value):
        self._max_retries = value

    @property
    def read_suff(self):
        return self._read_suff

    @read_suff.setter
    def read_suff(self, value):
        self._read_suff = value

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

    exec_args.pop("metadata", None)
    exec_args.update(kwargs)

    assembler = assembler_factory(cfg.asm.name, cfg.asm)
    max_attempts = assembler.max_retries
    for retry in range(max_attempts):
        try:
            assembler(**exec_args)
            break
        except Exception as ex:
            logger.info(f"Exception {ex} during execution of the attempt {retry + 1}/{max_attempts}")

    return get_job_outputs(exec_args)
