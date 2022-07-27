import tempfile
from pathlib import Path

from reads import reads_simulator


def test_pbmsim2_construct_exe_cmd(test_reads_root):
    cfg1 = {
        'name': 'pbsim2',
        'params': {
            'long': {
                'depth': 30,
            },
            'pattern': None
        }
    }

    cfg2 = {
        'name': 'pbsim2',
        'params': {
        }
    }

    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp_dir:
        vendor_dir = Path(tmp_dir) / 'vendor'
        vendor_dir.mkdir(exist_ok=True)
        simulator_root = vendor_dir / 'pbsim2'
        simulator_root.mkdir(exist_ok=True)

        dummy_ref = test_reads_root / 'dummy' / 'dummy.fa'
        save_dir = Path("/tmp/save")
        prefix = 'prefix'

        pbsim = reads_simulator.PbSim2(cfg1, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")} --depth 30 --prefix {prefix} {str(dummy_ref)}',
            f'mv {prefix}_0001.* {save_dir}'
        ]

        assert len(cmd) == 2
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]

        pbsim = reads_simulator.PbSim2(cfg2, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)
        assert len(cmd) == 2

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")}  --prefix {prefix} {str(dummy_ref)}',
            f'mv {prefix}_0001.* {save_dir}'
        ]
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]