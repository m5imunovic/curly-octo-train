from pathlib import Path

import pytest
from hydra import compose, initialize_config_dir


@pytest.fixture(scope="session")
def test_data_reads_fq(test_data_reads) -> list[Path]:
    return list((test_data_reads / "simulated").glob("*.fastq"))


@pytest.fixture(scope="session")
def test_la_jolla_cfg(test_cfg_root):
    test_asm_config_dir = test_cfg_root / "asm"
    with initialize_config_dir(version_base=None, config_dir=str(test_asm_config_dir), job_name="test_asm"):
        cfg = compose(config_name="test_la_jolla.yaml")
        return cfg
