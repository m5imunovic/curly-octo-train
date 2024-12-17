import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

from omegaconf import OmegaConf
from typeguard import typechecked

from reads.rsimulator import RSimulator
from utils import path_helpers as ph
from utils.io_utils import compose_cmd_params

OmegaConf.register_new_resolver("project_root", ph.project_root_append, replace=True)


logger = logging.getLogger(__name__)


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
            cmds.append(f"ln -s {profile_file} {profile_file.name}")
            cmds.append(f"ln -s {stats_file} {stats_file.name}")

        cmds.extend(
            [
                f"{self.simulator_exec} {option_params} {reference_params}",
            ]
        )

        if self.cfg.profile.path and self.cfg.params.long["sample-profile-id"]:
            cmds.append(f"rm {profile_file.name}")
            cmds.append(f"rm {stats_file.name}")

        return cmds

    @typechecked
    def run(self, genome: list[Path], output_path: Path, *args, **kwargs) -> bool:
        self._simulate_reads(genome, output_path)
        return True

    @typechecked
    def _simulate_reads(self, chr_seq_path: list[Path], output_path: Path):
        reference_list = "\n-->".join(str(p) for p in chr_seq_path)
        logger.info(f"Simulating reads from reference:\n-->{reference_list}")
        commands = self.construct_exec_cmd(chr_seq_path)
        with tempfile.TemporaryDirectory() as staging_dir:
            for cmd in commands:
                logger.info(f"RUN::simulate:: {cmd}")
                subprocess.run(cmd, shell=True, cwd=staging_dir)
            if not output_path.exists():
                output_path.mkdir(parents=True)
            # shutil.move would also include staging dir
            shutil.copytree(staging_dir, output_path, dirs_exist_ok=True, copy_function=shutil.move)
