import pytest
from hydra import compose, initialize_config_dir


@pytest.fixture
def expected_lja_dot():
    return {
        "number_of_nodes": 3502,
        "number_of_edges": 5116,
    }


@pytest.fixture(scope="session")
def test_db_graph_cfg(test_cfg_root):
    test_graph_dir = test_cfg_root / "graph"
    with initialize_config_dir(version_base=None, config_dir=str(test_graph_dir), job_name="test_graph"):
        cfg = compose(config_name="test_db_graph.yaml")
        return cfg
