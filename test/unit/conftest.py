import pytest
from omegaconf import OmegaConf

from utils.path_helpers import get_project_root, project_root_append

OmegaConf.register_new_resolver("project_root", project_root_append, replace=True)


@pytest.fixture(scope="session")
def test_data_root():
    return get_project_root() / "test" / "data"


@pytest.fixture(scope="session")
def test_reads_root(test_data_root):
    return test_data_root / "reads"


@pytest.fixture(scope="session")
def test_gfa_root(test_data_root):
    return test_data_root / "gfa"


@pytest.fixture(scope="session")
def test_graph_root(test_data_root):
    return test_data_root / "graph"


@pytest.fixture(scope="session")
def test_datasets_root(test_data_root):
    return test_data_root / "datasets"


@pytest.fixture(scope="session")
def test_data_reference(test_data_root):
    return test_data_root / "reference"


@pytest.fixture(scope="session")
def test_cfg_path(test_data_root):
    return test_data_root / "config" / "config.yaml"


@pytest.fixture(scope="session")
def test_cfg_root():
    return get_project_root() / "test" / "config"
