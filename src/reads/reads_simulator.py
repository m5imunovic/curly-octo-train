import multiprocessing as mp
import os
import shutil
import subprocess
from abc import abstractmethod
from pathlib import Path
from typing import List

import yaml
from omegaconf import DictConfig, OmegaConf
from typeguard import typechecked

from utils import path_helpers as ph
from utils.io_utils import compose_cmd_params, get_read_files

OmegaConf.register_new_resolver("project_root", ph.project_root_append, replace=True)


class RSimulator:
    def __init__(self, cfg: DictConfig, vendor_dir: Path):
        self.cfg = cfg
        self.name = self.cfg.name
        self.seed = self.cfg.params.long.seed or 42
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
        pass

    def pre_simulation_step(self, *args, **kwargs):
        pass

    def post_simulation_step(self, *args, **kwargs):
        pass

    @abstractmethod
    def run(self, ref_root: Path, simulated_data_root: Path, *args, **kwargs):
        pass


class PbSim2(RSimulator):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / "pbsim2"
        if not simulator_root.exists():
            try:
                print("SETUP::generate:: Download PbSim2")
                subprocess.run("git clone https://github.com/yukiteruono/pbsim2.git", shell=True, cwd=vendor_dir)
                subprocess.run("./configure", shell=True, cwd=simulator_root)
                subprocess.run("make", shell=True, cwd=simulator_root)
            except Exception as ex:
                print(f"SETUP::generate:: Error: {ex}")
                shutil.rmtree(simulator_root)
                raise ex
        else:
            print("SETUP::generate:: Use existing PbSim2")

        return simulator_root / "src" / "pbsim"

    def _install_from_script(self, script_path: Path):
        print(f"Skipping execution of {script_path}")
        raise NotImplementedError("Installation from script is not implemented!")

    @typechecked
    def _construct_exec_cmd(self, ref_path: Path, chr_save_path: Path, prefix: str) -> List[str]:
        assert "params" in self.cfg, "params must be specified in config"

        read_params = self._construct_read_params(ref_path)
        option_params = compose_cmd_params(self.cfg.params)
        prefix_param = f"--prefix {prefix}"
        read_params = self._prefix_read_params(read_params)

        cmds = []
        if self.cfg.profile_path and self.cfg.params.long["sample-profile-id"]:
            profile_id = self.cfg.params.long["sample-profile-id"]
            profile_file = Path(self.cfg.profile_path) / f"sample_profile_{profile_id}.fastq"
            stats_file = Path(self.cfg.profile_path) / f"sample_profile_{profile_id}.stats"
            assert profile_file.exists(), f"Profile file {profile_file} does not exist!"
            assert stats_file.exists(), f"Stats file {stats_file} does not exist!"
            cmds.append(f"ln -s {profile_file} {profile_file.name}")
            cmds.append(f"ln -s {stats_file} {stats_file.name}")

        cmds.extend(
            [
                f"{self.simulator_exec} {option_params} {prefix_param} {read_params}",
                f"rm {prefix}_0001.ref",
            ]
        )

        if self.cfg.profile_path and self.cfg.params.long["sample-profile-id"]:
            cmds.append(f"rm {profile_file.name}")
            cmds.append(f"rm {stats_file.name}")

        return cmds

    @typechecked
    def _construct_read_params(self, ref_path: Path):
        suffix = [".fasta", ".fa"] if "suffix" not in self.cfg else OmegaConf.to_container(self.cfg["suffix"])
        read_files = get_read_files(ref_path, suffix=suffix)
        read_params = " ".join(f"{str(read_file)}" for read_file in read_files)
        self._prefix_read_params(read_params)
        return read_params

    @typechecked
    def _prefix_read_params(self, read_params: str):
        return read_params

    @typechecked
    def run(self, ref_root: Path, simulated_species_path: Path, *args, **kwargs) -> bool:
        chr_path = ref_root / "chromosomes"
        assert chr_path.exists(), f"{chr_path} does not exist!"

        ref_suffix = list(self.cfg.suffix) or [".fasta", ".fa"]
        ref_fasta_files = get_read_files(chr_path, suffix=ref_suffix)
        simulation_data = []
        for ref_fasta_file in ref_fasta_files:
            simulated_fastq_file = simulated_species_path / f"{ref_fasta_file.stem}.fastq"
            simulation_data.append([ref_fasta_file, simulated_fastq_file, f"{ref_fasta_file.stem}"])

        print(f"SETUP::simulate:: Simulate {simulated_fastq_file} datasets from {ref_fasta_file}")
        # leave one processor free
        threads = self.cfg.threads or (os.cpu_count() - 1)
        with mp.Pool(threads) as pool:
            pool.starmap(self.simulate_reads_mp, simulation_data)

        return True

    @typechecked
    def simulate_reads_mp(self, chr_seq_path: Path, chr_save_path: Path, prefix: str):
        print(f"Simulating reads from reference:\n --> {chr_seq_path}")
        commands = self._construct_exec_cmd(chr_seq_path, chr_save_path, prefix)
        cwd_path = chr_save_path.parent
        if not cwd_path.exists():
            cwd_path.mkdir(parents=True)
        for cmd in commands:
            subprocess.run(cmd, shell=True, cwd=cwd_path)

    @typechecked
    def pre_simulation_step(self, simulated_species_path: Path, *args, **kwargs):
        if simulated_species_path.exists():
            if self.cfg.overwrite:
                print("PRE::simulate:: Removing existing simulation data")
                shutil.rmtree(simulated_species_path)
                # ensure that the directory exists
            else:
                err = f"PRE::simulate:: {simulated_species_path} already exists! Set overwrite to True to overwrite existing data!"
                raise FileExistsError(err)

        simulated_species_path.mkdir(parents=True)


class PbSim3(PbSim2):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / "pbsim3"
        if not simulator_root.exists():
            try:
                print("SETUP::generate:: Download PbSim3")
                subprocess.run("git clone https://github.com/yukiteruono/pbsim3.git", shell=True, cwd=vendor_dir)
                subprocess.run("./configure", shell=True, cwd=simulator_root)
                subprocess.run("make", shell=True, cwd=simulator_root)
            except Exception as ex:
                print("SETUP::generate:: Error: {ex}")
                shutil.rmtree(simulator_root)
                raise ex
        else:
            print("SETUP::generate:: Use existing PbSim3")

        return simulator_root / "src" / "pbsim"

    @typechecked
    def _prefix_read_params(self, read_params: str):
        return f"--genome {read_params}"


@typechecked
def simulator_factory(simulator: str, cfg: DictConfig) -> RSimulator:
    vendor_dir: Path = ph.get_vendor_path()
    if simulator == "pbsim2":
        return PbSim2(cfg=cfg, vendor_dir=vendor_dir)
    elif simulator == "pbsim3":
        return PbSim3(cfg=cfg, vendor_dir=vendor_dir)
    raise ValueError(f"Unknown simulator name {simulator}")


def run(cfg: DictConfig, **kwargs):
    output_path = Path(cfg.paths.reads_dir) / cfg.species_name["name"] / cfg.date_mm_dd / f"S{cfg.seed}"
    exec_args = {
        # Top level output path
        "simulated_species_path": output_path,
        # Path to the reference genome directory (can contain one or multiple fasta files)
        "ref_root": Path(cfg.paths.ref_dir) / cfg.species_name["name"],
        "experiment": cfg.experiment,
    }
    exec_args.update(kwargs)
    simulator = simulator_factory(simulator=cfg.reads.name, cfg=cfg.reads)

    simulator.pre_simulation_step(**exec_args)
    simulator.run(**exec_args)

    metadata_path = output_path.parent / "metadata"
    metadata_path.mkdir(parents=True, exist_ok=True)

    with open(metadata_path / f"S{cfg.seed}.yaml", "w") as f:
        OmegaConf.resolve(cfg)
        cfg_cont = OmegaConf.to_container(cfg)
        cfg_cont.pop("paths", None)
        yaml.dump(cfg_cont, f)
