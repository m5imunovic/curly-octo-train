"""
Adds additional features to the graph based on config
"""
from typing import List, Sequence, Set, Tuple

import networkx as nx
import numpy as np
from scipy import sparse
from typeguard import typechecked

from graph.construct_graph import DbGraphType


def supported_digraph_features() -> Set[str]:
    return {'ln', 'kc', 'in_degree', 'out_degree', 'pr_1', 'pr_2', 'pr_3', 'pr_4', 'pr_5'}


def add_digraph_features(g, features) -> Tuple[nx.DiGraph, List[str]]:
    available_features = ['kc', 'ln']
    if 'in_degree' in features:
        in_degrees = dict(g.in_degree())
        nx.set_node_attributes(g, in_degrees, 'in_degree')
        available_features.append('in_degree')
    if 'out_degree' in features:
        out_degrees = dict(g.out_degree())
        nx.set_node_attributes(g, out_degrees, 'out_degree')
        available_features.append('out_degree')

    for feature in features:
        # TODO: Does it makes sense to do this with torch instead and use GPU?
        if feature.startswith('pr_'):
            # get number of hops
            k = int(feature.split('_')[1])
            assert k > 0 and k < 6, f'Invalid number of hops {k}'
            # get adjacency matrix
            n = g.number_of_nodes()
            A = nx.linalg.adjacency_matrix(g)
            D = A.sum(axis=1)
            Dinv = 1./ (D+1e-9); Dinv[D<1e-9] = 0 # take care of nodes without outgoing edges
            Dinv = sparse.diags(np.squeeze(np.asarray(Dinv)), dtype=float) # D^-1 
            P = (Dinv @ A).T 
            pk_0 = np.ones(n) / n
            x = pk_0.copy()
            alpha = 0.95 
            for k_idx in range(k):
                x = alpha * P.dot(x) + (1-alpha) * pk_0
                attr = {node: x[idx] for idx, node in enumerate(g.nodes)}
                feature_name = f'pr_{k_idx+1}'
                nx.set_node_attributes(g, attr, feature_name)
                available_features.append(feature_name)
    
    return g, available_features


@typechecked
def add_features(g: DbGraphType, features: Sequence[str]) -> Tuple[DbGraphType, List[str]]:
    if isinstance(g, nx.MultiDiGraph):
        raise ValueError('Features are not yet supported for multigraphs')
    elif isinstance(g, nx.DiGraph):
        assert all(feat in supported_digraph_features() for feat in features)
        return add_digraph_features(g, features)
        # assert supported features
    else:
        raise ValueError(f'Unknown graph type {type(g)}')
