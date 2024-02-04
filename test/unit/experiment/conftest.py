from pathlib import Path

import pytest
from hydra import compose, initialize_config_dir

from experiment.scenario_schema import Scenario

this_dir = Path(__file__).parent


@pytest.fixture(scope="session")
def test_data_experiment(test_data_root) -> Path:
    return test_data_root / "experiment"


@pytest.fixture(scope="session")
def test_scenario(test_cfg_root) -> Scenario:
    test_scenario_path = test_cfg_root / "experiment" / "scenarios" / "test_scenario.json"
    with open(test_scenario_path) as f:
        scenario = Scenario.from_json(f.read())

    return scenario


@pytest.fixture(scope="session")
def test_experiment_cfg(test_cfg_root):
    test_experiment_config_dir = test_cfg_root / "experiment"
    with initialize_config_dir(version_base=None, config_dir=str(test_experiment_config_dir), job_name="test_exp"):
        cfg = compose(config_name="test_experiment.yaml")
        return cfg
