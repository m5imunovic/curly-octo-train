import os
from pathlib import Path

import pytest
from omegaconf import OmegaConf

from utils.path_helpers import get_project_root, project_root_append

OmegaConf.register_new_resolver("project_root", project_root_append, replace=True)


@pytest.fixture(scope="session", autouse=True)
def set_env():
    os.environ["TEST_DATA_ROOT"] = str(get_project_root() / "test" / "data")
    os.environ["TEST_PROFILE_ROOT"] = str(get_project_root() / "test" / "data" / "pf")


@pytest.fixture(scope="session")
def test_data_root():
    return get_project_root() / "test" / "data"


@pytest.fixture(scope="session")
def test_reads_root(test_data_root):
    return test_data_root / "reads"


@pytest.fixture(scope="session")
def test_graph_root(test_data_root):
    return test_data_root / "graph"


@pytest.fixture(scope="session")
def test_dot_root(test_graph_root):
    return test_graph_root / "dot"


@pytest.fixture(scope="session")
def test_datasets_root(test_data_root):
    return test_data_root / "datasets"


@pytest.fixture(scope="session")
def test_data_reference(test_data_root):
    return test_data_root / "references"


@pytest.fixture(scope="session")
def test_cfg_path(test_data_root):
    return test_data_root / "config" / "config.yaml"


@pytest.fixture(scope="session")
def test_cfg_root():
    return get_project_root() / "test" / "config"


@pytest.fixture(scope="session")
def test_species_genome(test_data_reference) -> list[Path]:
    chr_path = test_data_reference / "test_species" / "S_0001" / "chromosomes"
    return list(chr_path.glob("*.fasta"))


@pytest.fixture(scope="session")
def test_species_reads_root(test_data_reference) -> list[Path]:
    chr_path = test_data_reference / "test_species" / "S_0001" / "reads"
    return chr_path


@pytest.fixture(scope="session")
def test_data_reads(test_data_root) -> Path:
    return test_data_root / "reads"


@pytest.fixture(scope="session")
def test_data_assemblies(test_data_root) -> Path:
    return test_data_root / "assemblies"
