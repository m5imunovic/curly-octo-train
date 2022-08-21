from typing import Dict, List, Optional, Union

import hydra
import networkx as nx
from omegaconf import DictConfig
from torch_geometric.data import Data
from torch_geometric.utils.convert import from_networkx
from typeguard import typechecked

from graph.construct_graph import construct_graph, DbGraphType
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
        y_idx = pyg.node_attr.shape[1] - 1
        pyg.y = pyg.node_attr[:, y_idx]
        pyg.node_attr = pyg.node_attr[:, :y_idx]
    return pyg


@typechecked
def convert_to_pyg_graph(g: DbGraphType, has_mult_info: bool = False) -> Data:
    group_attrs = ['kc', 'ln']
    if has_mult_info:
        # supervised data available
        group_attrs.append('mult')
    if isinstance(g, nx.MultiDiGraph):
        return convert_to_pyg_multigraph(g, group_edge_attrs=group_attrs, has_mult_info=has_mult_info)
    if isinstance(g, nx.DiGraph):
        return convert_to_pyg_digraph(g, group_node_attrs=group_attrs, has_mult_info=has_mult_info)

    raise KeyError(f'Graph type {type(g)} not supported')


def run(cfg: DictConfig):
    g = construct_graph(cfg)
    print(f'Number of edges {g.number_of_edges()}')
    print(f'Number of nodes {g.number_of_nodes()}')
    mult_info = parse_mult_info(cfg.mult_info_path)
    g = add_mult_info_features(g, mult_info)
    pyg = convert_to_pyg_graph(g, has_mult_info=True)


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()