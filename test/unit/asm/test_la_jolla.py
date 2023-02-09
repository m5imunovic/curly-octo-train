import tempfile
from pathlib import Path

import omegaconf

from asm import la_jolla


def test_la_jolla_construct_exe_cmd(test_reads_root):
    cfg1 = omegaconf.DictConfig({
        'name': 'LJA',
        'params': {
            'long': {
                'threads': 8
            },
            'short': {
                'k': 501
            }, 
            'append':  '--compress'
        },
        'exec': 'jumboDBG',
        'full_asm': 'full_asm',
    })

    with tempfile.TemporaryDirectory() as tmp_dir:
        vendor_dir = Path(tmp_dir) / 'vendor'
        vendor_dir.mkdir(exist_ok=True)
        assembler_root = vendor_dir / 'LJA'
        assembler_root.mkdir(exist_ok=True)

        lja = la_jolla.LaJolla(cfg1, vendor_dir)
        reads_path = test_reads_root / 'chr1'
        output_path = Path(tmp_dir) / 'output'
        reads = [reads_path / '0_0001.fastq', reads_path / '1_0001.fastq']
        cmd = lja._construct_exec_cmd(reads_files=reads, ref_path=reads_path, output_path=output_path)

        # call jumboDBG
        # call lja
        # call align_and_print twice
        assert len(cmd) == 4

        assert str(assembler_root / 'bin/jumboDBG') in cmd[0]
        assert '-k 501' in cmd[0]
        assert '--threads 8' in cmd[0]
        assert '--compress' in cmd[0]
        assert '--reference' in cmd[0]
        assert str(reads_path / '0_0001.fastq') in cmd[0]
        assert str(reads_path / '1_0001.fastq') in cmd[0]
        assert str(assembler_root / 'bin/lja') in cmd[1]
        assert '--compress' not in cmd[1]
        assert '--reference' not in cmd[1]