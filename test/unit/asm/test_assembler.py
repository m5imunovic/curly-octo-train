from tty import CFLAG
from unittest import mock

from omegaconf import OmegaConf

from asm.assembler import run as run_assembler
from asm.la_jolla import LaJolla


def assembler_cfg(overwrite: bool):
    return OmegaConf.create({
        'name': 'LJA',
        'overwrite': overwrite,
        'params': None,
        'experiment': 'test_experiment',
        'species': 'test_species'
    })


@mock.patch.object(LaJolla, 'run')
@mock.patch.object(LaJolla, '_install')
@mock.patch('shutil.rmtree', return_value=True)
def test_assembler_overwrite_data(mock_rmtree, mock_install, mock_run, tmp_path):
    cfg = assembler_cfg(overwrite=True)

    kwargs = {
        'reads_path': tmp_path,
        'out_path': tmp_path
    }

    mock_install.return_value = tmp_path
    mock_run.return_value = True
    run_assembler(cfg, **kwargs)

    assert mock_rmtree.call_count == 1
    assert mock_run.call_count == 1

    cfg = assembler_cfg(overwrite=False)
    run_assembler(cfg, **kwargs)
    assert mock_rmtree.call_count == 1
    assert mock_run.call_count == 2