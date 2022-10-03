# TODO:
# test that the order (i.e. adjacency operations are OK) is preserved
# test with real data
import pytest

import networkx as nx
from omegaconf import OmegaConf

from graph.construct_graph import construct_graph
from graph.construct_features import add_features


@pytest.fixture
def simple_digraph() -> nx.DiGraph:
    # nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (4, 7), (5, 8), (6, 8)]
    return nx.DiGraph(incoming_graph_data=edges)


@pytest.fixture
def simple_multidigraph() -> nx.MultiDiGraph:
    # nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    edges = [(0, 1), (0,1), (0, 2), (1, 3), (1, 4), (1, 4), (2, 5), (2, 6), (3, 7), (3, 7), (3,7), (4, 7), (5, 8), (6, 8)]
    return nx.MultiDiGraph(incoming_graph_data=edges)


def test_add_features_in_degree_out_degree(simple_digraph):
    g, _ = add_features(simple_digraph, ['in_degree', 'out_degree'])
    expected_in_degrees = dict(g.in_degree())
    expected_out_degrees = dict(g.out_degree())
    for node_with_data in g.nodes(data=True):
        nid, attrs = node_with_data
        expected_in_degree = expected_in_degrees[nid]
        expected_out_degree = expected_out_degrees[nid]
        assert expected_in_degree == attrs['in_degree']
        assert expected_out_degree == attrs['out_degree']


def test_add_features_in_degree_out_degree_multigraph(simple_multidigraph):
    g, _ = add_features(simple_multidigraph, ['in_degree', 'out_degree'])
    expected_in_degrees = dict(g.in_degree())
    expected_out_degrees = dict(g.out_degree())
    for node_with_data in g.nodes(data=True):
        nid, attrs = node_with_data
        expected_in_degree = expected_in_degrees[nid]
        expected_out_degree = expected_out_degrees[nid]
        assert expected_in_degree == attrs['in_degree']
        assert expected_out_degree == attrs['out_degree']


def test_add_features_pr(simple_digraph):
    g, _ = add_features(simple_digraph, ['pr_1'])
    expected_pr1 = {0: 0.005555, 1: 0.058333, 2: 0.058333, 3: 0.058333, 4: 0.058333, 5: 0.058333, 6: 0.058333, 7: 0.216666, 8: 0.216666}
    for node_with_data in g.nodes(data=True):
        nid, attrs = node_with_data
        assert pytest.approx(expected_pr1[nid], rel=1e-3) == attrs['pr_1']

    g, _ = add_features(simple_digraph, ['pr_4'])
    for node_with_data in g.nodes(data=True):
        _, attrs = node_with_data
        assert 'pr_1' in attrs
        assert 'pr_2' in attrs
        assert 'pr_3' in attrs
        assert 'pr_4' in attrs


def test_add_features_pr_multidigraph(simple_multidigraph):
    g, _ = add_features(simple_multidigraph, ['pr_1'])
    expected_pr1 = {0: 0.005555, 1: 0.075925, 2: 0.040740, 3: 0.040740, 4: 0.075925, 5: 0.058333, 6: 0.058333, 7: 0.216666, 8: 0.216666}
    for node_with_data in g.nodes(data=True):
        nid, attrs = node_with_data
        assert pytest.approx(expected_pr1[nid], rel=1e-3) == attrs['pr_1']

    g, _ = add_features(simple_multidigraph, ['pr_4'])
    for node_with_data in g.nodes(data=True):
        _, attrs = node_with_data
        assert 'pr_1' in attrs
        assert 'pr_2' in attrs
        assert 'pr_3' in attrs
        assert 'pr_4' in attrs


#@pytest.mark.parametrize('graph_type', ['digraph', 'multigraph'])
def test_add_features_lja_graph(test_gfa_root): #, graph_type):
    graph_type = 'digraph'
    cfg = OmegaConf.create({
        'graph_type': graph_type, 
        'gfa_path': test_gfa_root / 'lja_graph.gfa',
        'k': 501
    })

    g = construct_graph(cfg)
    g, _ = add_features(g, features=['ln', 'kc', 'pr_4', 'in_degree', 'out_degree'])
    for node_with_data in g.nodes(data=True):
        _, attrs = node_with_data
        assert 'ln' in attrs
        assert 'kc' in attrs
        assert 'pr_1' in attrs
        assert 'pr_2' in attrs
        assert 'pr_3' in attrs
        assert 'pr_4' in attrs
        assert 'in_degree' in attrs
        assert 'out_degree' in attrs