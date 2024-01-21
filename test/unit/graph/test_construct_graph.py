import pytest
from omegaconf import OmegaConf

from graph.construct_graph import construct_graphs


@pytest.mark.parametrize("graph_type", ["digraph", "multidigraph"])
def test_construct_example1_graph(test_gfa_root, graph_type):
    cfg = OmegaConf.create({"graph_type": graph_type, "gfa_path": test_gfa_root / "example1.gfa", "k": 10})

    with pytest.raises(ValueError) as e:
        construct_graphs(cfg.gfa_path, cfg.k)


def test_construct_lja_graph(test_gfa_root, expected_lja_dict):
    cfg = OmegaConf.create({"gfa_path": test_gfa_root / "lja_graph.gfa", "k": 501})

    g, labels = construct_graphs(cfg.gfa_path, cfg.k)
    for g_type in ["multidigraph"]:  # , "digraph"]:
        assert g_type in g
        assert g_type in labels
        assert g[g_type].number_of_nodes() == expected_lja_dict[g_type]["number_of_nodes"]
        assert g[g_type].number_of_edges() == expected_lja_dict[g_type]["number_of_edges"]
