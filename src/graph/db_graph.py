import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

import networkx as nx
import torch
from omegaconf import DictConfig
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from typeguard import typechecked

from asm.mult_info_parser import parse_mult_info
from graph.construct_features import FeatureDict, add_features
from graph.construct_graph import DbGraphType, construct_graph

logger = logging.getLogger(__name__)


def add_mult_info_features(g: DbGraphType, mult_info: Dict[str, int]) -> DbGraphType:
    # We use a small trick here. As we name this feature 'y' it automatically gets assigned to the y
    # attribute of the exported pytorch geometric graph and can then be used for learning supervision
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]["y"] = mult_info[key]
    else:
        raise TypeError(f"Graph type {type(g)} not supported")
    return g


def convert_to_pyg_multigraph(g: DbGraphType, group_attrs: Optional[FeatureDict] = None) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_attrs["node"], group_edge_attrs=group_attrs["edge"])
    return pyg


@typechecked
def convert_to_pyg_graph(g: DbGraphType, features: Optional[FeatureDict] = None) -> Data:
    if isinstance(g, nx.MultiDiGraph):
        return convert_to_pyg_multigraph(g, group_attrs=features)

    raise TypeError(f"Graph type {type(g)} not supported")


def export_ids(g: DbGraphType, out_path: str) -> None:
    id_map = {}
    if isinstance(g, nx.MultiDiGraph):
        for idx, edge in enumerate(g.edges(keys=True)):
            id_map[idx] = edge[2]
    else:
        raise TypeError(f"Graph type {type(g)} not supported")

    with open(out_path, "w") as handle:
        json.dump(id_map, handle, indent=4)


def process_graph(idx: int, assembly_path: Path, cfg: DictConfig, out_raw_path: Path, out_debug_paths: Path):
    # Process raw data
    mult_info_path = assembly_path / "mult.info"
    gfa_path = assembly_path / "graph.gfa"
    logger.info(f"Processing {gfa_path}")

    result = {"raw_files": defaultdict(tuple)}
    g, labels = construct_graph(gfa_path=gfa_path, k=cfg.k)

    g, features = add_features(g, cfg.features)
    logger.info(f"{idx}: Number of edges {g.number_of_edges()}")
    logger.info(f"{idx}: Number of nodes {g.number_of_nodes()}")
    mult_info = parse_mult_info(mult_info_path)
    g = add_mult_info_features(g, mult_info)
    # TODO: only ever output these only for the test data
    if cfg.debug:
        nx.write_graphml(g, out_debug_paths / f"{idx}.graphml")

    with open(out_debug_paths / f"{idx}.rcmap", "w") as f:
        json.dump(labels, f, indent=4)

    pyg = convert_to_pyg_graph(g, features)
    pyg_filename = f"{idx}.pt"
    logger.info(f"Saving {out_raw_path / pyg_filename }")
    torch.save(pyg, out_raw_path / pyg_filename)
    result["raw_files"] = (pyg_filename, assembly_path)

    export_ids(g, out_debug_paths / f"{idx}.idmap")

    return result


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        "assembly_path": None,
        "output_path": None,
        "idx": 0,
    }

    exec_args.update(kwargs)

    # Iterate over experiment directory to find individual assemblies

    output_path = exec_args["output_path"]
    # TODO:
    # We should add collect and clean step which would save the files in a way we want
    out_raw_path = output_path / "raw"
    out_debug_path = output_path / "debug"

    if not out_raw_path.exists():
        out_raw_path.mkdir(parents=True)
    if not out_debug_path.exists():
        out_debug_path.mkdir(parents=True)

    idx = exec_args["idx"]
    assembly_path = exec_args["assembly_path"]
    processed_results = process_graph(idx, assembly_path, cfg, out_raw_path, out_debug_path)

    with open(output_path / "raw.csv", "w") as f:
        pyg_file, origin_asm_dir = processed_results["raw_files"]
        f.write(f"{pyg_file, origin_asm_dir}\n")
