import pytest

from utils.path_helpers import get_project_root

@pytest.fixture()
def test_reads_root():
    return get_project_root() / 'test' / 'data' / 'reads'


@pytest.fixture()
def test_gfa_root():
    return get_project_root() / 'test' / 'data' / 'gfa'


