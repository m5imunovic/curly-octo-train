import subprocess
import shutil
from pathlib import Path

from omegaconf import OmegaConf
from typeguard import typechecked

from asm import assembler
from utils.io_utils import compose_cmd_params, get_read_files


class LaJolla(assembler.Assembler):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        assembler_root = vendor_dir / 'LJA'
        if not assembler_root.exists():
            try:
                print(f'SETUP::generate:: Download La Jolla Assembler')
                subprocess.run('git clone https://github.com/AntonBankevich/LJA.git', shell=True, cwd=str(vendor_dir))
                subprocess.run('git checkout development', shell=True, cwd=str(assembler_root))
                subprocess.run('cmake .', shell=True, cwd=str(assembler_root))
                subprocess.run('make -j 8', shell=True, cwd=str(assembler_root))
            except Exception as ex:
                print(f'SETUP::generate:: Error: {ex}')
                shutil.rmtree(assembler_root)
                raise ex

        return assembler_root / 'bin'

    @typechecked
    def _construct_exec_cmd(self, reads_path: Path, output_path: Path) -> list[str]:
        assert 'params' in self.cfg, "params must be specified in config"

        suffix = ['.fastq', '.fq'] if 'suffix' not in self.cfg else OmegaConf.to_container(self.cfg['suffix'])
        read_files = get_read_files(reads_path, suffix=suffix)
        reads_cmd_params = ' '.join([f'--reads {str(read_file)}' for read_file in read_files])
        out_cmd_param = f'-o {str(output_path)}'
        option_params = compose_cmd_params(self.cfg['params'])

        exec = 'lja' if 'exec' not in self.cfg else self.cfg['exec']
        asm_executable = self.assembler_root / exec
        return [f'{asm_executable} {option_params} {reads_cmd_params} {out_cmd_param}']

    @typechecked
    def run(self, reads_path: Path, out_path: Path, *args, **kwargs):
        commands = self._construct_exec_cmd(reads_path, out_path)
        for cmd in commands:
            print(f'RUN::assembler::\n{cmd}')
            subprocess.run(cmd, shell=True)

    @typechecked
    def pre_assembly_step(self, out_path: Path, *args, **kwargs):
        if self.cfg.overwrite:
            print(f'PRE::assembler:: Remove existing assembly files')
            if out_path.exists():
                shutil.rmtree(out_path)
                out_path.mkdir(parents=True, exist_ok=True)

