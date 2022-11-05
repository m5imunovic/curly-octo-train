import shutil
from collections import defaultdict
from itertools import chain
from typing import Dict, Optional

import hydra
import networkx as nx
import torch
from omegaconf import DictConfig
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from typeguard import typechecked

import utils.path_helpers as ph
from graph.construct_features import add_features, FeatureDict
from graph.construct_graph import construct_graphs, DbGraphType
from graph.mult_info_parser import parse_mult_info


def add_mult_info_features(g: DbGraphType, mult_info: Dict[str, int]) -> DbGraphType:
    # We use a small trick here. As we name this feature 'y' it automatically gets assigned to the y
    # attribute of the exported pytorch geometric graph and can then be used for learning supervision
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]['y'] = mult_info[key]
    else:
        nx.set_node_attributes(g, mult_info, 'y')
    return g


def convert_to_pyg_multigraph(g: DbGraphType, group_attrs: Optional[FeatureDict] = None) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_attrs['node'], group_edge_attrs=group_attrs['edge'])
    return pyg


def convert_to_pyg_digraph(g: DbGraphType, group_attrs: Optional[FeatureDict] = None) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_attrs['node'], group_edge_attrs=group_attrs['edge'])
    return pyg


@typechecked
def convert_to_pyg_graph(g: DbGraphType, features: Optional[FeatureDict] = None) -> Data:
    if isinstance(g, nx.MultiDiGraph):
        return convert_to_pyg_multigraph(g, group_attrs=features)
    if isinstance(g, nx.DiGraph):
        return convert_to_pyg_digraph(g, group_attrs=features)

    raise KeyError(f'Graph type {type(g)} not supported')


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        'assemblies_path': cfg.paths.assemblies_dir / cfg.asm.experiment,
        'out_path': cfg.paths.datasets_dir / cfg.graph.experiment
    }

    exec_args.update(kwargs)

    # Iterate over experiment directory to find individual assemblies
    assemblies_path = exec_args['assemblies_path']
    chr_dirs = list(filter(lambda x: x.is_dir(), assemblies_path.iterdir()))
    graph_dirs = [list(filter(lambda x: x.is_dir() and x.stem.isdigit(), chr_dir.iterdir())) for chr_dir in chr_dirs]
    graph_dirs = chain(*graph_dirs)

    out_path = exec_args['out_path']
    out_raw_paths = {}
    out_processed_paths = {}
    out_debug_paths = {}
    for g_type in ['digraph', 'multidigraph']:
        out_raw_path = out_path / g_type / 'raw'
        out_processed_path = out_path / g_type / 'processed'
        out_debug_path = out_path / g_type / 'debug'

        if not out_raw_path.exists():
            out_raw_path.mkdir(parents=True)
        if not out_processed_path.exists():
            out_processed_path.mkdir(parents=True)
        if not out_debug_path.exists():
            out_debug_path.mkdir(parents=True)
        out_raw_paths[g_type] = out_raw_path
        out_processed_paths[g_type] = out_processed_path
        out_debug_paths[g_type] = out_debug_path

    processed_files = defaultdict(list)
    raw_files = defaultdict(list)

    # TODO: parallelize (checko out joblib package)
    for idx, graph_dir in enumerate(graph_dirs):
        # Process raw data
        cfg.graph.mult_info_path = graph_dir / 'mult.info'
        cfg.graph.gfa_path = graph_dir / 'graph.gfa'

        all_graphs = construct_graphs(cfg.graph)
        for g_type, g in all_graphs.items():
            g, features = add_features(g, cfg.graph.features)
            print(f'Number of edges {g.number_of_edges()}')
            print(f'Number of nodes {g.number_of_nodes()}')
            mult_info = parse_mult_info(cfg.graph.mult_info_path)
            g = add_mult_info_features(g, mult_info)
            if cfg.graph.debug:
                nx.write_gml(g, out_debug_paths[g_type] / f'{idx}.gml')
                
            pyg = convert_to_pyg_graph(g, features)
            pyg_filename = f'{idx}.pt'
            torch.save(pyg, out_processed_paths[g_type] / pyg_filename)
            processed_files[g_type].append((pyg_filename, graph_dir))

            # Save raw data
            raw_filename = f'{idx}'
            shutil.make_archive(out_raw_paths[g_type] / raw_filename, 'zip', graph_dir)
            raw_files[g_type].append((raw_filename, graph_dir))


    for g_type, p_files in processed_files.items():
        with open(out_path / g_type / 'processed.csv', 'w') as f:
            for file in p_files:
                f.write(f'{file[0]},{file[1]}\n')

    for g_type, r_files in raw_files.items():
        with open(out_path / g_type / 'raw.csv', 'w') as f:
            for file in r_files:
                f.write(f'{file[0]},{file[1]}\n')


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()