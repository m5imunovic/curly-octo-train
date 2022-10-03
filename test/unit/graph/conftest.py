import pytest

import networkx as nx


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
