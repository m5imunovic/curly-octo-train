import pytest

import networkx as nx
from omegaconf import OmegaConf
from graph.construct_graph import construct_graph


expected_result = {
    'digraph': {
        'instance': nx.DiGraph,
        'number_of_nodes': 40,
        'number_of_edges': 24,
    },
    'multigraph': {
        'instance': nx.MultiDiGraph,
        'number_of_nodes': 40,
        'number_of_edges': 24,
    }
}


@pytest.fixture
def expected(request):
    graph_type = request.node.funcargs['graph_type']
    return expected_result[graph_type]


@pytest.mark.parametrize('graph_type', ['digraph', pytest.param('multigraph', marks=pytest.mark.xfail)])
def test_construct_graph(test_gfa_root, graph_type, expected):
    cfg = OmegaConf.create({
        'graph_type': graph_type, 
        'gfa_path': test_gfa_root / 'example1.gfa',
        'k': 10
    })

    g = construct_graph(cfg)
    assert isinstance(g, expected['instance'])
    assert g.number_of_nodes() == expected['number_of_nodes']
    assert g.number_of_edges() == expected['number_of_edges']