import subprocess
import shutil
from pathlib import Path
from typing import List

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
    def _construct_exec_cmd(self, reads_path: Path, output_path: Path) -> List[str]:
        assert 'params' in self.cfg, "params must be specified in config"

        suffix = ['.fastq', '.fq'] if 'suffix' not in self.cfg else OmegaConf.to_container(self.cfg['suffix'])
        read_files = get_read_files(reads_path, suffix=suffix)
        reads_cmd_params = ' '.join([f'--reads {str(read_file)}' for read_file in read_files])
        out_cmd_param = f'-o {str(output_path)}'
        option_params = compose_cmd_params(self.cfg['params'])

        exec = 'lja' if 'exec' not in self.cfg else self.cfg['exec']
        asm_executable = self.assembler_root / exec
        cmds = [f'{asm_executable} {option_params} {reads_cmd_params} {out_cmd_param}']
        return cmds


    @typechecked
    def run(self, reads_path: Path, out_path: Path, *args, **kwargs):
        for path in reads_path.iterdir():
            if path.is_dir():
                stem = path.stem
                commands = self._construct_exec_cmd(reads_path / stem, out_path / stem)
                cwd_path = out_path / stem
                if not cwd_path.exists():
                    cwd_path.mkdir(parents=True)
                for cmd in commands:
                    print(f'RUN::assembler::\n{cmd}')
                    subprocess.run(cmd, shell=True, cwd=cwd_path)

    @typechecked
    def pre_assembly_step(self, out_path: Path, *args, **kwargs):
        if self.cfg.overwrite:
            print(f'PRE::assembler:: Remove existing assembly files')
            if out_path.exists():
                shutil.rmtree(out_path)
                out_path.mkdir(parents=True, exist_ok=True)


    def post_assembly_step(self, out_path: Path, *args, **kwargs):
        # Clean up the files
        if 'keep' in self.cfg and self.cfg['keep'] is not None:
            keep_files = OmegaConf.to_container(self.cfg['keep'])
            for path in out_path.glob('**/*'):
                if not path.is_dir():
                    if path.name not in keep_files:
                        path.unlink()
        # Create files summary and store it to the csv file in raw directory
        kept_files = []
        for path in out_path.glob('**/*'):
            if not path.is_dir():
                kept_files.append(str(path.relative_to(out_path)))
        
        raw_dir = out_path / 'raw'
        if not raw_dir.exists():
            raw_dir.mkdir(parents=True)
        
        with open(raw_dir / 'files.csv', 'w') as f:
            f.write('\n'.join(kept_files))


