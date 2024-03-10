from pathlib import Path

from omegaconf import OmegaConf

from asm.mult_info_parser import parse_mult_info
from graph.construct_features import add_mult_info_features
from graph.db_graph import convert_to_pyg_graph, run
from graph.dot_parser import custom_parse_dot


def test_convert_to_pyg_graph_dot(test_dot_root, expected_lja_dot):
    cfg = OmegaConf.create({"graph_path": test_dot_root / "example1.dot", "k": 501})

    g = custom_parse_dot(cfg.graph_path, cfg.k)

    features = {"node": None, "edge": ["kc", "ln"]}
    mult_info_path = test_dot_root / "example1_mult.info"
    mult_info = parse_mult_info(mult_info_path)
    g = add_mult_info_features(g, mult_info=mult_info)
    pyg = convert_to_pyg_graph(g, features)
    expected_number_of_nodes = expected_lja_dot["number_of_nodes"]
    expected_number_of_edges = expected_lja_dot["number_of_edges"]
    assert pyg.num_nodes == expected_number_of_nodes
    assert pyg.edge_index.shape[1] == expected_number_of_edges
    assert pyg.y.shape[0] == expected_number_of_edges


def test_graph_produces_expected_outputs(test_db_graph_cfg, test_data_assemblies, tmpdir):
    exec_args = {
        "assembly_path": test_data_assemblies,
        "output_path": Path(tmpdir),
        "idx": 0,
    }

    run(test_db_graph_cfg, **exec_args)

    idx = exec_args["idx"]
    assert (Path(tmpdir) / "raw" / f"{idx}.pt").exists()
    assert (Path(tmpdir) / "debug" / f"{idx}.idmap").exists()
