import pytest

from utils.path_helpers import get_project_root


@pytest.fixture(scope='session')
def test_data_root():
    return get_project_root() / 'test' / 'data'


@pytest.fixture(scope='session')
def test_reads_root(test_data_root):
    return test_data_root / 'reads'


@pytest.fixture(scope='session')
def test_gfa_root(test_data_root):
    return test_data_root / 'gfa'


@pytest.fixture(scope='session')
def test_graph_root(test_data_root):
    return test_data_root / 'graph'


@pytest.fixture(scope='session')
def test_datasets_root(test_data_root):
    return test_data_root / 'datasets'