import pytest

from omegaconf import OmegaConf

from graph.construct_graph import construct_graph
from graph.db_graph import convert_to_pyg_graph


@pytest.mark.parametrize('graph_type', ['digraph', 'multidigraph'])
def test_convert_to_pyg_graph(test_gfa_root, graph_type, expected_lja):
    cfg = OmegaConf.create({
        'graph_type': graph_type, 
        'gfa_path': test_gfa_root / 'lja_graph.gfa',
        'k': 501
    })

    g, labels = construct_graph(cfg)

    if graph_type == 'digraph':
        features = {'node': ['kc', 'ln'], 'edge': None}

    else:
        features = {'node': None, 'edge': ['kc', 'ln']}

    expected_number_of_nodes = expected_lja['number_of_nodes']
    expected_number_of_edges = expected_lja['number_of_edges']


    pyg = convert_to_pyg_graph(g, features)

    assert pyg.num_nodes == expected_number_of_nodes
    assert pyg.edge_index.shape[1] == expected_number_of_edges