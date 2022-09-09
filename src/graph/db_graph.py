import shutil
from itertools import chain
from typing import Dict, List, Optional, Union

import hydra
import networkx as nx
import torch
from omegaconf import DictConfig, OmegaConf
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from typeguard import typechecked

import utils.path_helpers as ph
from graph.construct_features import add_features
from graph.construct_graph import construct_graph, DbGraphType
from graph.dbg_dataset import DBGDataset
from graph.mult_info_parser import parse_mult_info


FeatureList = Union[List[str], all]


def add_mult_info_features(g: DbGraphType, mult_info: Dict[str, int]) -> DbGraphType:
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]['mult'] = mult_info[key]
    else:
        nx.set_node_attributes(g, mult_info, 'mult')
    return g


def convert_to_pyg_multigraph(g: DbGraphType, group_edge_attrs: Optional[FeatureList] = None, has_mult_info: bool = True) -> Data:
    pyg = from_networkx(g, group_edge_attrs=group_edge_attrs)
    if has_mult_info:
        y_idx = pyg.edge_attr.shape[1] - 1
        pyg.y = pyg.edge_attr[:, y_idx]
        pyg.edge_attr = pyg.edge_attr[:, :y_idx]
    return pyg


def convert_to_pyg_digraph(g: DbGraphType, group_node_attrs: Optional[FeatureList] = None, has_mult_info: bool = True) -> Data:
    pyg = from_networkx(g, group_node_attrs=group_node_attrs)
    if has_mult_info:
        y_idx = pyg.x.shape[1] - 1
        pyg.y = pyg.x[:, y_idx]
        pyg.x = pyg.x[:, :y_idx]
    return pyg


@typechecked
def convert_to_pyg_graph(g: DbGraphType, group_attrs: List,  has_mult_info: bool = False) -> Data:
    if has_mult_info:
        # supervised data available
        group_attrs.append('mult')
    if isinstance(g, nx.MultiDiGraph):
        return convert_to_pyg_multigraph(g, group_edge_attrs=group_attrs, has_mult_info=has_mult_info)
    if isinstance(g, nx.DiGraph):
        return convert_to_pyg_digraph(g, group_node_attrs=group_attrs, has_mult_info=has_mult_info)

    raise KeyError(f'Graph type {type(g)} not supported')


def run(cfg: DictConfig, **kwargs):
    exec_args = {
        'assemblies_path': ph.get_assemblies_path() / cfg.experiment,
        'out_path': ph.get_datasets_path() / cfg.experiment
    }

    exec_args.update(kwargs)

    # Iterate over experiment directory to find individual assemblies
    assemblies_path = exec_args['assemblies_path']
    chr_dirs = list(filter(lambda x: x.is_dir(), assemblies_path.iterdir()))
    graph_dirs = [list(filter(lambda x: x.is_dir() and x.stem.isdigit(), chr_dir.iterdir())) for chr_dir in chr_dirs]
    graph_dirs = chain(*graph_dirs)

    out_path = exec_args['out_path']
    out_raw_path = out_path / 'raw'
    out_processed_path = out_path / 'processed'
    if not out_raw_path.exists():
        out_raw_path.mkdir(parents=True)
    if not out_processed_path.exists():
        out_processed_path.mkdir(parents=True)

    processed_files = []
    raw_files = []

    # TODO: parallelize (checko out joblib package)
    for idx, graph_dir in enumerate(graph_dirs):
        # Process raw data
        cfg.mult_info_path = graph_dir / 'mult.info'
        cfg.gfa_path = graph_dir / 'graph.gfa'

        g = construct_graph(cfg)
        g, features = add_features(g, cfg.features)
        print(f'Number of edges {g.number_of_edges()}')
        print(f'Number of nodes {g.number_of_nodes()}')
        mult_info = parse_mult_info(cfg.mult_info_path)
        g = add_mult_info_features(g, mult_info)
        pyg = convert_to_pyg_graph(g, features, has_mult_info=True)
        pyg_filename = f'{idx}.pt'
        torch.save(pyg, out_processed_path / pyg_filename)
        processed_files.append((pyg_filename, graph_dir))

        # Save raw data
        raw_filename = f'{idx}'
        shutil.make_archive(out_raw_path / raw_filename, 'zip', graph_dir)
        raw_files.append((raw_filename, graph_dir))

    with open(out_path / 'processed.csv', 'w') as f:
        for file in processed_files:
            f.write(f'{file[0]},{file[1]}\n')

    with open(out_path/ 'raw.csv', 'w') as f:
        for file in raw_files:
            f.write(f'{file[0]},{file[1]}\n')


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()