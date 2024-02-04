from pathlib import Path

import pytest
from hydra import compose, initialize_config_dir

this_dir = Path(__file__).parent


@pytest.fixture(scope="session")
def test_data_reads(test_data_root) -> Path:
    return test_data_root / "reads"


@pytest.fixture(scope="session")
def test_data_profile(test_data_root) -> Path:
    return test_data_root / "pf"


@pytest.fixture(scope="session")
def test_pbsim3_pf_cfg(test_cfg_root, test_data_profile):
    test_species_config_dir = test_cfg_root / "reads"
    overrides = [f"profile.path={test_data_profile}", "params.long.seed=2"]
    with initialize_config_dir(version_base=None, config_dir=str(test_species_config_dir), job_name="test_seq"):
        cfg = compose(config_name="test_pbsim3_pf.yaml", overrides=overrides)
        return cfg
