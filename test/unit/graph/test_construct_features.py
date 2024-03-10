# TODO:
# test that the order (i.e. adjacency operations are OK) is preserved
# test with real data
import networkx as nx
import pytest
from omegaconf import OmegaConf

from graph.construct_features import add_features
from graph.dot_parser import custom_parse_dot


def test_add_features_lja_graph_gfa(test_dot_root):
    cfg = OmegaConf.create({"graph_path": test_dot_root / "example1.dot", "k": 501})

    g = custom_parse_dot(cfg.graph_path, cfg.k)
    with pytest.raises(NotImplementedError):
        g, _ = add_features(g, features=["ln", "kc"])
    for node_with_data in g.nodes(data=True):
        _, attrs = node_with_data
        assert len(attrs) == 0
    for edge_with_data in g.edges(data=True):
        _, _, attrs = edge_with_data
        assert "ln" in attrs
        assert "kc" in attrs
