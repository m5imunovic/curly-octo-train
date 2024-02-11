from pathlib import Path

from omegaconf import OmegaConf

from graph.construct_graph import construct_graph
from graph.db_graph import convert_to_pyg_graph, run


def test_convert_to_pyg_graph(test_gfa_root, expected_lja):
    cfg = OmegaConf.create({"gfa_path": test_gfa_root / "lja_graph.gfa", "k": 501})

    g, _ = construct_graph(cfg.gfa_path, cfg.k)

    features = {"node": None, "edge": ["kc", "ln"]}

    expected_number_of_nodes = expected_lja["number_of_nodes"]
    expected_number_of_edges = expected_lja["number_of_edges"]

    pyg = convert_to_pyg_graph(g, features)

    assert pyg.num_nodes == expected_number_of_nodes
    assert pyg.edge_index.shape[1] == expected_number_of_edges


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
    assert (Path(tmpdir) / "debug" / f"{idx}.rcmap").exists()
