import pytest
from omegaconf import OmegaConf

from graph.construct_graph import construct_graph


def test_construct_example1_graph_error(test_gfa_root):
    cfg = OmegaConf.create({"gfa_path": test_gfa_root / "example1.gfa", "k": 10})

    with pytest.raises(ValueError) as e:
        construct_graph(cfg.gfa_path, cfg.k)


def test_construct_lja_graph(test_gfa_root, expected_lja):
    cfg = OmegaConf.create({"gfa_path": test_gfa_root / "lja_graph.gfa", "k": 501})

    g, _ = construct_graph(cfg.gfa_path, cfg.k)
    assert g.number_of_nodes() == expected_lja["number_of_nodes"]
    assert g.number_of_edges() == expected_lja["number_of_edges"]
