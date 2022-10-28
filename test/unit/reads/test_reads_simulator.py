import os
import tempfile
from pathlib import Path
from unittest import mock

from omegaconf import OmegaConf

from reads import reads_simulator
from reads.reads_simulator import PbSim2
from reads.reads_simulator import run as run_reads_simulator


def reads_simulator_cfg(overwrite: bool):
    return OmegaConf.create({
            'name': 'pbsim2',
            'overwrite': overwrite,
            'params': None,
            'species': 'test_species',
            'request': {
                'chr1': 2
            }
    })


@mock.patch('random.randint', return_value=22)
def test_pbmsim2_construct_exe_cmd(mock_randint, test_reads_root):
    cfg1 = OmegaConf.create({
        'name': 'pbsim2',
        'params': {
            'long': {
                'depth': 30,
            },
            'pattern': None
        }
    })

    cfg2 = OmegaConf.create({
        'name': 'pbsim2',
        'params': {
        }
    })

    with tempfile.TemporaryDirectory() as tmp_dir:
        vendor_dir = Path(tmp_dir) / 'vendor'
        vendor_dir.mkdir(exist_ok=True)
        simulator_root = vendor_dir / 'pbsim2'
        simulator_root.mkdir(exist_ok=True)

        dummy_ref = test_reads_root / 'dummy' / 'dummy.fa'
        save_dir = Path("/tmp/save")
        prefix = '1'
        pbsim = reads_simulator.PbSim2(cfg1, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)
        seed = mock_randint.return_value + int(prefix)

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")} --depth 30 --prefix {prefix} {str(dummy_ref)} --seed {seed}',
            f'rm {prefix}_0001.ref'
        ]

        assert len(cmd) == 2
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]

        pbsim = reads_simulator.PbSim2(cfg2, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)
        assert len(cmd) == 2

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")}  --prefix {prefix} {str(dummy_ref)} --seed {seed}',
            f'rm {prefix}_0001.ref'
        ]
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]


@mock.patch.object(PbSim2, 'run')
@mock.patch('shutil.rmtree', return_value=True)
def test_reads_simulator_overwrite_data(mock_rmtree, mock_run, tmp_path):
    cfg = reads_simulator_cfg(overwrite=True)
    kwargs = {
        'ref_root': tmp_path / cfg['species'],
        'simulated_data_root': tmp_path,
        'chr_request': dict(cfg['request'])
    }
    mock_run.return_value = True

    os.mkdir(str(tmp_path / cfg['species']))
    run_reads_simulator(cfg, **kwargs)
    assert mock_rmtree.call_count == 1
    assert mock_run.call_count == 1

    cfg = reads_simulator_cfg(overwrite=False)
    run_reads_simulator(cfg, **kwargs)
    assert mock_rmtree.call_count == 1
    assert mock_run.call_count == 2