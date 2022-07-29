import tempfile
from pathlib import Path

from asm import la_jolla


def test_la_jolla_construct_exe_cmd(test_reads_root):
    cfg1 = {
        'name': 'LJA',
        'params': {
            'long': {
                'threads': 8
            },
            'short': {
                'k': 501
            }, 
            'append':  '--compress'
        }
    }

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
        vendor_dir = Path(tmp_dir) / 'vendor'
        vendor_dir.mkdir(exist_ok=True)
        assembler_root = vendor_dir / 'LJA'
        assembler_root.mkdir(exist_ok=True)

        lja = la_jolla.LaJolla(cfg1, vendor_dir)
        reads_path = test_reads_root / 'chr1'
        output_path = Path(tmp_dir) / 'output'
        cmd = lja._construct_exec_cmd(reads_path, output_path)

        assert len(cmd) == 1
        print(cmd)

        assert str(assembler_root / 'bin/lja') in cmd[0]
        assert '-k 501 --threads 8 --compress' in cmd[0]
        assert str(reads_path / '0_0001.fastq') in cmd[0]
        assert str(reads_path / '1_0001.fastq') in cmd[0]