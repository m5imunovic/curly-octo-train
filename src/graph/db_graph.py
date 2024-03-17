import json
import logging
from pathlib import Path
from typing import Optional

import networkx as nx
import torch
from omegaconf import DictConfig
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from typeguard import typechecked

from asm.mult_info_parser import parse_mult_info
from graph.construct_features import FeatureDict, add_features, add_mult_info_features
from graph.dot_parser import custom_parse_dot
from utils.io_utils import get_job_outputs

logger = logging.getLogger(__name__)
DbGraphType = nx.MultiDiGraph


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


def process_graph(idx: int, assembly_path: Path, cfg: DictConfig, output_path: Path):
    out_raw_path = output_path / "raw"
    out_debug_path = output_path / "debug"

    if not out_raw_path.exists():
        out_raw_path.mkdir(parents=True)
    if not out_debug_path.exists():
        out_debug_path.mkdir(parents=True)

    # Process raw data
    graph_path = assembly_path / "graph.dot"
    logger.info(f"Processing {graph_path}...")
    g = custom_parse_dot(graph_path, k=cfg.k)
    logger.info(f"{idx}: Number of edges {g.number_of_edges()}")
    logger.info(f"{idx}: Number of nodes {g.number_of_nodes()}")

    logger.info("Loading multiplicity info...")
    mult_info_path = assembly_path / "mult.info"
    mult_info = parse_mult_info(mult_info_path)
    g, features = add_features(g, cfg.features)
    g = add_mult_info_features(g, mult_info)
    logger.info("Converting to torch format...")
    pyg = convert_to_pyg_graph(g, features)
    pyg_filename = f"{idx}.pt"
    logger.info(f"Saving {out_raw_path / pyg_filename }")
    torch.save(pyg, out_raw_path / pyg_filename)
    if cfg.debug:
        nx.write_graphml(g, out_debug_path / f"{idx}.graphml")

    export_ids(g, out_debug_path / f"{idx}.idmap")


def run(cfg: DictConfig, **kwargs) -> dict:
    exec_args = {
        "assembly_path": None,
        "output_path": None,
    }

    exec_args.update(kwargs)

    output_path = exec_args["output_path"]
    assembly_path = exec_args["assembly_path"]
    idx = 0
    process_graph(idx, assembly_path, cfg, output_path)

    return get_job_outputs(exec_args)
