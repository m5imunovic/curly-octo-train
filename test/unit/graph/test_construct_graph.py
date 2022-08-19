import pytest

import networkx as nx
from omegaconf import OmegaConf
from graph.construct_graph import construct_graph


expected_lja_graph = {
    'digraph': {
        'instance': nx.DiGraph,
        'number_of_nodes': 160,
        'number_of_edges': 66,
    },
    'multigraph': {
        'instance': nx.MultiDiGraph,
        'number_of_nodes': 188,
        'number_of_edges': 160,
    }
}


@pytest.fixture
def expected_lja(request):
    graph_type = request.node.funcargs['graph_type']
    return expected_lja_graph[graph_type]


@pytest.mark.parametrize('graph_type', ['digraph', 'multigraph'])
def test_construct_example1_graph(test_gfa_root, graph_type):
    cfg = OmegaConf.create({
        'graph_type': graph_type, 
        'gfa_path': test_gfa_root / 'example1.gfa',
        'k': 10
    })

    with pytest.raises(ValueError):
        construct_graph(cfg)


@pytest.mark.parametrize('graph_type', ['digraph', 'multigraph'])
def test_construct_lja_graph(test_gfa_root, graph_type, expected_lja):
    cfg = OmegaConf.create({
        'graph_type': graph_type, 
        'gfa_path': test_gfa_root / 'lja_graph.gfa',
        'k': 501
    })

    g = construct_graph(cfg)
    assert isinstance(g, expected_lja['instance'])
    assert g.number_of_nodes() == expected_lja['number_of_nodes']
    assert g.number_of_edges() == expected_lja['number_of_edges']