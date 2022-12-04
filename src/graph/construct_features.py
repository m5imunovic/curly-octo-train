"""
Adds additional features to the graph based on config
"""
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

import networkx as nx
import numpy as np
from scipy import sparse
from typeguard import typechecked

from graph.construct_graph import DbGraphType


FeatureDict = Union[Dict[str, Optional[List]], all]


def supported_digraph_features() -> Set[str]:
    return {'ln', 'kc'}


def supported_multidigraph_features() -> Set[str]:
    return {'ln', 'kc'}


@typechecked
def add_multidigraph_features(g: nx.MultiDiGraph, features: Sequence[str]) -> Tuple[nx.MultiDiGraph, FeatureDict]:
    available_edge_features = ['kc', 'ln']
    available_node_features = None

    return g, {'edge': available_edge_features, 'node': available_node_features}


@typechecked
def add_digraph_features(g: nx.DiGraph, features: Sequence[str]) -> Tuple[nx.DiGraph, FeatureDict]:
    available_node_features = ['kc', 'ln']
    
    return g, {'node': available_node_features, 'edge': None}


@typechecked
def add_features(g: DbGraphType, features: Sequence[str]) -> Tuple[DbGraphType, FeatureDict]:
    if isinstance(g, nx.MultiDiGraph):
        assert all(feat in supported_multidigraph_features() for feat in features)
        return add_multidigraph_features(g, features)
    elif isinstance(g, nx.DiGraph):
        assert all(feat in supported_digraph_features() for feat in features)
        return add_digraph_features(g, features)
    else:
        raise ValueError(f'Unknown graph type {type(g)}')
