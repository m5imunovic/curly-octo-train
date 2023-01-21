import subprocess
import shutil
from collections import defaultdict
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
    def _construct_exec_cmd(self, reads_files: List[Path], ref_path: Path, output_path: Path) -> List[str]:
        assert 'params' in self.cfg, "params must be specified in config"

        reads_cmd_params = ' '.join([f'--reads {str(read_file)}' for read_file in reads_files])
        self.cfg.params.long.reference = ref_path
        out_cmd_param = f'-o {str(output_path)}'
        option_params = compose_cmd_params(self.cfg['params'])


        exec = 'lja' if 'exec' not in self.cfg else self.cfg['exec']
        asm_executable = self.assembler_root / exec
        cmds = [f'{asm_executable} {option_params} {reads_cmd_params} {out_cmd_param}']

        if self.cfg['full_asm']:
            lja_executable = self.assembler_root / 'lja'
            full_asm_subdir = self.cfg['full_asm']
            params = OmegaConf.to_container(self.cfg['params'])
            params['long'].pop('reference', None)
            params['append'] = params['append'].replace('--compress', '')
            params['short'].update({'K': 5001})

            option_params = compose_cmd_params(params)
            out_cmd_param = f'-o {str(output_path / full_asm_subdir)}'

            cmds.append(f'{lja_executable} {option_params} {reads_cmd_params} {out_cmd_param}')
        return cmds


    @typechecked
    def run(self, ref_path: Path, reads_path: Path, out_path: Path, suffix, *args, **kwargs):
        # TODO: check if num ref files smaller than reads files
        ref_paths = get_read_files(ref_path, suffix=suffix['reference'])

        fastq_files = get_read_files(reads_path, suffix=suffix['reads'], regex=self.cfg.dir_filter)
        fastq_groups = defaultdict(list)
        for fastq_file in fastq_files:
            group = fastq_file.parent.relative_to(reads_path)
            fastq_groups[group].append(fastq_file)
 

        fasta_stems = {ref_path.stem: ref_path for ref_path in ref_paths}
        for group, fastq_files in fastq_groups.items():
            # TODO: add regex stem matching
            fastq_stems = {fastq_file.stem[:-5] : fastq_file for fastq_file in fastq_files}
            for fastq_stem in fastq_stems:
                commands = self._construct_exec_cmd(
                    [fastq_stems[fastq_stem]], fasta_stems[fastq_stem], out_path / group / fastq_stem
                )
            
                cwd_path = out_path / group / fastq_stem
                if not cwd_path.exists():
                    cwd_path.mkdir(parents=True)
                if self.cfg['full_asm']:
                    (cwd_path / self.cfg['full_asm']).mkdir(parents=True)
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


