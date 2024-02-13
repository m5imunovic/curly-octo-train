import pytest
from hydra import compose, initialize_config_dir


@pytest.fixture(scope="session")
def test_pbsim3_pf_cfg(test_cfg_root):
    test_reads_config_dir = test_cfg_root / "reads"
    overrides = ["params.long.seed=2"]
    with initialize_config_dir(version_base=None, config_dir=str(test_reads_config_dir), job_name="test_seq"):
        cfg = compose(config_name="test_pbsim3_pf.yaml", overrides=overrides)
        return cfg
