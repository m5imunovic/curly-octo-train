from typing import Dict

import hydra
import networkx as nx
from omegaconf import DictConfig

from graph.construct_graph import construct_graph, DbGraphType
from graph.mult_info_parser import parse_mult_info


def add_mult_info_features(g: DbGraphType, mult_info: Dict[str, int]) -> DbGraphType:
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]['mult'] = mult_info[key]
    else:
        nx.set_node_attributes(g, mult_info, 'mult')
    return g


def run(cfg: DictConfig):
    g = construct_graph(cfg)
    print(f'Number of edges {g.number_of_edges()}')
    print(f'Number of nodes {g.number_of_nodes()}')
    mult_info = parse_mult_info(cfg.mult_info_path)
    g = add_mult_info_features(g, mult_info)


@hydra.main(version_base=None, config_path='../../config/graph', config_name='db_graph')
def main(cfg: DictConfig):
    run(cfg)


if __name__ == '__main__':
    main()