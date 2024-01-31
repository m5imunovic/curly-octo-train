from pathlib import Path

import pytest
from hydra import compose, initialize_config_dir

this_dir = Path(__file__).parent


@pytest.fixture(scope="session")
def test_data_reference(test_data_root):
    return test_data_root / "reference"


@pytest.fixture(scope="session")
def test_species_cfg(test_cfg_root):
    test_species_config_dir = test_cfg_root / "reference" / "species"
    with initialize_config_dir(version_base=None, config_dir=str(test_species_config_dir), job_name="test_ref"):
        cfg = compose(config_name="test_species.yaml")
        return cfg
