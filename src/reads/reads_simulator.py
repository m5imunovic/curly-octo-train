import multiprocessing as mp
import os
import shutil
import subprocess
from abc import abstractmethod
from pathlib import Path

import hydra
from omegaconf import OmegaConf
from typeguard import typechecked

from utils import path_helpers as ph
from utils.io_utils import compose_cmd_params, get_read_files


OmegaConf.register_new_resolver('project_root', ph.project_root_append)


class RSimulator:
    def __init__(self, cfg, vendor_dir: Path):
        self.cfg = cfg
        if 'install_script' in cfg and cfg['install_script'] is not None:
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
    def run(self, ref_root: Path, simulated_data_root: Path, chr_request: dict, *args, **kwargs):
        pass


class PbSim2(RSimulator):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        simulator_root = vendor_dir / 'pbsim2'
        if not simulator_root.exists():
            try:
                print(f'SETUP::generate:: Download PbSim2')
                subprocess.run('git clone https://github.com/yukiteruono/pbsim2.git', shell=True, cwd=vendor_dir)
                subprocess.run(f'./configure', shell=True, cwd=simulator_root)
                subprocess.run(f'make', shell=True, cwd=simulator_root)
            except Exception as ex:
                print(f'SETUP::generate:: Error: {ex}')
                shutil.rmtree(simulator_root)
                raise ex
        else:
            print(f'SETUP::generate:: Use existing PbSim2')

        return simulator_root / 'src' / 'pbsim'

    def _install_from_script(self, script_path: Path):
        print(f'Skipping execution of {script_path}')
        raise NotImplementedError('Installation from script is not implemented!')

    @typechecked
    def _construct_exec_cmd(self, ref_path: Path, chr_save_path: Path, prefix: str) -> list[str]:
        assert 'params' in self.cfg, "params must be specified in config"

        suffix = ['.fasta', '.fa'] if 'suffix' not in self.cfg else OmegaConf.to_container(self.cfg['suffix'])
        read_files = get_read_files(ref_path, suffix=suffix)
        read_params = ' '.join(f'{str(read_file)}' for read_file in read_files)
        option_params = compose_cmd_params(self.cfg['params'])
        prefix_param = f'--prefix {prefix}'

        return [
            f'{self.simulator_exec} {option_params} {prefix_param} {read_params}',
            f'mv {prefix}_0001.* {chr_save_path}'
        ]

    @typechecked
    def run(self, ref_root: Path, simulated_data_root: Path, chr_request: dict, *args, **kwargs):
        chr_path = ref_root / 'chromosomes'
        assert chr_path.exists(), f'{chr_path} does not exist!'
        assert simulated_data_root.exists(), f'{simulated_data_root} does not exist!'
        simulation_data = []
        for chrN, n_need in chr_request.items():
            species_name = ref_root.stem
            chr_raw_path = simulated_data_root / species_name / f'{chrN}'
            if not chr_raw_path.exists():
                chr_raw_path.mkdir(parents=True)
                n_have = 0
            else:
                n_have = len(get_read_files(chr_raw_path, suffix=['.fastq', '.fq']))
            if n_need <= n_have:
                continue
            else:
                n_diff = n_need - n_have
                print(f'SETUP::simulate:: Simulate {n_diff} datasets for {chrN}')
                # Simulate reads for chrN n_diff times

                # TODO: use get_read_files() and check if the respective chromosome reference exists
                chr_seq_path = chr_path / f'{chrN}.fasta'
                for i in range(n_diff):
                    idx = n_have + i
                    chr_save_path = chr_raw_path
                    simulation_data.append((chr_seq_path, chr_save_path, str(idx), '/'.join([str(i+1), str(n_diff)])))

            # leave one processor free
            with mp.Pool(os.cpu_count() - 1) as pool:
                pool.starmap(self.simulate_reads_mp, simulation_data)

    @typechecked
    def simulate_reads_mp(self, chr_seq_path: Path, chr_save_path: Path, prefix: str, i: str):
        print(f'Request {i}: Simulating reads from referece:\n --> {chr_seq_path}')
        commands = self._construct_exec_cmd(chr_seq_path, chr_save_path, prefix)
        for cmd in commands:
            subprocess.run(cmd, shell=True)


@typechecked
def simulator_factory(simulator: str, cfg: dict) -> RSimulator:
    vendor_dir: Path = ph.get_vendor_path()
    if simulator == 'pbsim2':
        return PbSim2(cfg=cfg, vendor_dir=vendor_dir)
    raise ValueError(f"Unknown simulator name {simulator}")


def run(cfg, **kwargs):
    read_simulator = dict(cfg)
    simulator = simulator_factory(simulator=read_simulator['name'], cfg=read_simulator)

    kwargs = {
        'ref_root': ph.get_ref_path() / read_simulator['species'],
        'simulated_data_root': ph.get_simulated_data_path(),
        'chr_request': dict(cfg['request'])
    }
    simulator.run(**kwargs)


@hydra.main(version_base=None, config_path='../../config/reads', config_name='pbsim2')
def main(cfg):
    print("Running read simulator step...")
    run(cfg)


if __name__ == '__main__':
    main()