from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig
from typeguard import typechecked

READ_FILE = "sim_0001.fastq"


class RSimulator:
    def __init__(self, cfg: DictConfig, vendor_dir: Path):
        self.cfg = cfg
        self.name = self.cfg.name
        if "install_script" in cfg and cfg.install_script is not None:
            assert cfg.exec_root is not None
            self._install_from_script(Path(cfg.script_path))
            self.simulator_exec = Path(cfg.exec_root)
        else:
            self.simulator_exec = self._install(vendor_dir)

    @typechecked
    @abstractmethod
    def _install(self, vendor_dir: Path):
        pass

    @typechecked
    def _install_from_script(self, script_path: Path):
        raise NotImplementedError(f"Ignoring {script_path} commands. Installation from script is not implemented!")

    def pre_simulation_step(self, *args, **kwargs):
        pass

    def post_simulation_step(self, *args, **kwargs):
        pass

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        self.pre_simulation_step(*args, **kwargs)
        self.run(*args, **kwargs)
        self.post_simulation_step(*args, **kwargs)
