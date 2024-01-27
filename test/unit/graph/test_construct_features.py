# TODO:
# test that the order (i.e. adjacency operations are OK) is preserved
# test with real data
import networkx as nx
import pytest
from omegaconf import OmegaConf

from graph.construct_features import add_features
from graph.construct_graph import construct_graph


@pytest.fixture
def simple_multidigraph() -> nx.MultiDiGraph:
    # nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    edges = [
        (0, 1),
        (0, 1),
        (0, 2),
        (1, 3),
        (1, 4),
        (1, 4),
        (2, 5),
        (2, 6),
        (3, 7),
        (3, 7),
        (3, 7),
        (4, 7),
        (5, 8),
        (6, 8),
    ]
    return nx.MultiDiGraph(incoming_graph_data=edges)


def test_add_features_lja_graph(test_gfa_root):
    cfg = OmegaConf.create({"gfa_path": test_gfa_root / "lja_graph.gfa", "k": 501})

    g, _ = construct_graph(cfg.gfa_path, cfg.k)
    with pytest.raises(NotImplementedError):
        g, _ = add_features(g, features=["ln", "kc"])
    for node_with_data in g.nodes(data=True):
        _, attrs = node_with_data
        assert len(attrs) == 0
    for edge_with_data in g.edges(data=True):
        _, _, attrs = edge_with_data
        assert "ln" in attrs
        assert "kc" in attrs
