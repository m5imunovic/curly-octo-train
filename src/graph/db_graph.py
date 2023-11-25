import json
import multiprocessing as mp
import os
import shutil
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
from graph.construct_graph import DbGraphType, construct_graphs
from utils.io_utils import get_read_files


def add_mult_info_features(g: DbGraphType, mult_info: Dict[str, int]) -> DbGraphType:
    # We use a small trick here. As we name this feature 'y' it automatically gets assigned to the y
    # attribute of the exported pytorch geometric graph and can then be used for learning supervision
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]["y"] = mult_info[key]
    else:
        nx.set_node_attributes(g, mult_info, "y")
    return g


def convert_to_pyg_multigraph(g: DbGraphType, group_attrs: Optional[FeatureDict] = None) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_attrs["node"], group_edge_attrs=group_attrs["edge"])
    return pyg


def convert_to_pyg_digraph(g: DbGraphType, group_attrs: Optional[FeatureDict] = None) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_attrs["node"], group_edge_attrs=group_attrs["edge"])
    return pyg


@typechecked
def convert_to_pyg_graph(g: DbGraphType, features: Optional[FeatureDict] = None) -> Data:
    if isinstance(g, nx.MultiDiGraph):
        return convert_to_pyg_multigraph(g, group_attrs=features)
    if isinstance(g, nx.DiGraph):
        return convert_to_pyg_digraph(g, group_attrs=features)

    raise KeyError(f"Graph type {type(g)} not supported")


def export_ids(g: DbGraphType, out_path: str) -> None:
    id_map = {}
    if isinstance(g, nx.MultiDiGraph):
        for idx, edge in enumerate(g.edges(keys=True)):
            id_map[idx] = edge[2]
    else:
        for idx, node in enumerate(g.nodes()):
            id_map[idx] = node

    with open(out_path, "w") as handle:
        json.dump(id_map, handle, indent=4)


def process_graph_mp(
    idx: int, graph_dir: Path, cfg: DictConfig, out_processed_paths: Path, out_raw_paths: Path, out_debug_paths: Path
):
    # Process raw data
    mult_info_path = graph_dir / "mult.info"
    gfa_path = graph_dir / "graph.gfa"
    print(f"Processing {gfa_path}")

    result = {"processed_files": defaultdict(tuple), "raw_files": defaultdict(tuple)}
    all_graphs, all_labels = construct_graphs(gfa_path=gfa_path, k=cfg.graph.k)

    for g_type, g in all_graphs.items():
        g, features = add_features(g, cfg.graph.features)
        print(f"{idx}: Number of edges {g.number_of_edges()}")
        print(f"{idx}: Number of nodes {g.number_of_nodes()}")
        mult_info = parse_mult_info(mult_info_path)
        g = add_mult_info_features(g, mult_info)
        if cfg.graph.debug:
            nx.write_graphml(g, out_debug_paths[g_type] / f"{idx}.graphml")

        with open(out_debug_paths[g_type] / f"{idx}.rcmap", "w") as f:
            json.dump(all_labels[g_type], f, indent=4)

        pyg = convert_to_pyg_graph(g, features)
        pyg_filename = f"{idx}.pt"
        print(f"Saving {out_processed_paths[g_type] / pyg_filename }")
        result["processed_files"][g_type] = (pyg_filename, graph_dir)

        export_ids(g, out_debug_paths[g_type] / f"{idx}.idmap")
        torch.save(pyg, out_processed_paths[g_type] / pyg_filename)

        # Save raw data
        raw_filename = f"{idx}"
        shutil.make_archive(out_raw_paths[g_type] / raw_filename, "zip", graph_dir)
        result["raw_files"][g_type] = (raw_filename, graph_dir)

    return result


def run(cfg: DictConfig, **kwargs):

    exec_args = {
        "assemblies_path": cfg.paths.assemblies_dir / cfg.species_name.name,
        "out_path": cfg.paths.datasets_dir / f"{cfg.species_name.name}_{cfg.date_mm_dd}" / cfg.graph.set.name,
    }

    exec_args.update(kwargs)

    # Iterate over experiment directory to find individual assemblies
    assemblies_path = exec_args["assemblies_path"]
    graph_dirs = [p.parent for p in get_read_files(assemblies_path, suffix=[".gfa"], regex=cfg.graph.set.dir_filter)]

    out_path = exec_args["out_path"]
    out_raw_paths = {}
    out_processed_paths = {}
    out_debug_paths = {}
    for g_type in ["digraph", "multidigraph"]:
        out_raw_path = out_path / g_type / "raw"
        out_processed_path = out_path / g_type / "processed"
        out_debug_path = out_path / g_type / "debug"

        if not out_raw_path.exists():
            out_raw_path.mkdir(parents=True)
        if not out_processed_path.exists():
            out_processed_path.mkdir(parents=True)
        if not out_debug_path.exists():
            out_debug_path.mkdir(parents=True)
        out_raw_paths[g_type] = out_raw_path
        out_processed_paths[g_type] = out_processed_path
        out_debug_paths[g_type] = out_debug_path

    # Prepare multiprocessing data
    mp_data = []
    threads = cfg.graph.threads or (os.cpu_count() - 1)
    for idx, graph_dir in enumerate(graph_dirs):
        mp_data.append([idx, graph_dir, cfg, out_processed_paths, out_raw_paths, out_debug_paths])

    print(f"Processing {len(mp_data)} graphs with {threads} threads")
    #if len(mp_data) == 1:
    #    processed_results = [process_graph_mp(*mp_data[0])]
    with mp.Pool(threads) as pool:
        processed_results = pool.starmap(process_graph_mp, mp_data)
    

    # Merge results
    processed_files = defaultdict(list)
    raw_files = defaultdict(list)
    for result in processed_results:
        for g_type, p_files in result["processed_files"].items():
            processed_files[g_type].append(p_files)
        for g_type, r_files in result["raw_files"].items():
            raw_files[g_type].append(r_files)

    for g_type, p_files in processed_files.items():
        with open(out_path / g_type / "processed.csv", "w") as f:
            for file in p_files:
                f.write(f"{file[0]},{file[1]}\n")

    for g_type, r_files in raw_files.items():
        with open(out_path / g_type / "raw.csv", "w") as f:
            for file in r_files:
                f.write(f"{file[0]},{file[1]}\n")
