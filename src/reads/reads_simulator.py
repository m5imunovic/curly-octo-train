import logging
import shutil
import subprocess
import tempfile
from abc import abstractmethod
from pathlib import Path

from omegaconf import DictConfig, OmegaConf
from typeguard import typechecked

from utils import path_helpers as ph
from utils.io_utils import compose_cmd_params

OmegaConf.register_new_resolver("project_root", ph.project_root_append, replace=True)


logger = logging.getLogger(__name__)


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


class PbSim3(RSimulator):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / "pbsim3"
        if not simulator_root.exists():
            try:
                logger.info("SETUP::generate:: Download PbSim3")
                subprocess.run("git clone https://github.com/yukiteruono/pbsim3.git", shell=True, cwd=vendor_dir)
                # TODO: hardcode checkout commit
                subprocess.run("./configure", shell=True, cwd=simulator_root)
                subprocess.run("make", shell=True, cwd=simulator_root)
            except Exception as ex:
                logger.error(f"SETUP::generate:: Error: {ex}")
                shutil.rmtree(simulator_root)
                raise ex
        else:
            logger.info("SETUP::generate:: Use existing PbSim3")

        return simulator_root / "src" / "pbsim"

    @typechecked
    def _construct_reference_params(self, genome: list[Path]):
        return " ".join(f"--genome {chr}" for chr in genome)

    @staticmethod
    def construct_sample_path(profile_root: Path | None, profile_id: str | None, suffix: str) -> Path | None:
        if profile_root and profile_id:
            return (profile_root / f"sample_profile_{profile_id}").with_suffix(suffix)
        return None

    @property
    def sample_profile_path(self) -> Path | None:
        return self.construct_sample_path(
            Path(self.cfg.profile.path), self.cfg.params.long["sample-profile-id"], ".fastq"
        )

    @property
    def sample_stats_path(self) -> Path | None:
        return self.construct_sample_path(
            Path(self.cfg.profile.path), self.cfg.params.long["sample-profile-id"], ".stats"
        )

    @typechecked
    def construct_exec_cmd(self, genome: list[Path]) -> list[str]:
        assert "params" in self.cfg, "params must be specified in config"

        reference_params = self._construct_reference_params(genome)
        option_params = compose_cmd_params(self.cfg.params)

        # TODO: add named command for prettier output
        cmds = []
        if self.cfg.profile.path and self.cfg.params.long["sample-profile-id"]:
            profile_file = self.sample_profile_path
            stats_file = self.sample_stats_path
            cmds.append(f"ln {profile_file} {profile_file.name}")
            cmds.append(f"ln {stats_file} {stats_file.name}")

        prefix = self.cfg.params.long.prefix
        cmds.extend(
            [
                f"{self.simulator_exec} {option_params} {reference_params}",
                f"rm {prefix}_0001.maf",
                f"rm {prefix}_0001.ref",
            ]
        )

        if self.cfg.profile.path and self.cfg.params.long["sample-profile-id"]:
            cmds.append(f"rm {profile_file.name}")
            cmds.append(f"rm {stats_file.name}")

        return cmds

    @typechecked
    def run(self, genome: list[Path], simulated_reads_path: Path, *args, **kwargs) -> bool:
        self._simulate_reads(genome, simulated_reads_path)
        return True

    @typechecked
    def _simulate_reads(self, chr_seq_path: list[Path], simulated_reads_path: Path):
        reference_list = "\n-->".join(str(p) for p in chr_seq_path)
        logger.info(f"Simulating reads from reference:\n-->{reference_list}")
        commands = self.construct_exec_cmd(chr_seq_path)
        with tempfile.TemporaryDirectory() as staging_dir:
            for cmd in commands:
                logger.info(f"RUN::simulate:: {cmd}")
                subprocess.run(cmd, shell=True, cwd=staging_dir)
            if not simulated_reads_path.exists():
                simulated_reads_path.mkdir(parents=True)
            # shutil.move would also include staging dir
            shutil.copytree(staging_dir, simulated_reads_path, dirs_exist_ok=True, copy_function=shutil.move)


@typechecked
def simulator_factory(simulator: str, cfg: DictConfig) -> RSimulator:
    vendor_dir: Path = ph.get_vendor_path()
    if simulator == "pbsim3":
        return PbSim3(cfg=cfg, vendor_dir=vendor_dir)
    raise ValueError(f"Unknown simulator name {simulator}")
