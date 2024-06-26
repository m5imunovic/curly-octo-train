import pytest
from hydra import compose, initialize_config_dir


@pytest.fixture(scope="session")
def test_species_cfg(test_cfg_root):
    test_species_config_dir = test_cfg_root / "reference" / "species"
    with initialize_config_dir(version_base=None, config_dir=str(test_species_config_dir), job_name="test_ref"):
        cfg = compose(config_name="test_species.yaml")
        return cfg
