"""Adds additional features to the graph based on config."""
from typing import Dict, List, Optional, Sequence, Set, Tuple, Union

import networkx as nx
from typeguard import typechecked

from graph.construct_graph import DbGraphType

FeatureDict = Union[Dict[str, Optional[List]], all]


def supported_multidigraph_features() -> Set[str]:
    return {"ln", "kc"}


@typechecked
def add_multidigraph_features(
    g: nx.MultiDiGraph, features: Sequence[str] | None
) -> Tuple[nx.MultiDiGraph, FeatureDict]:
    available_edge_features = ["kc", "ln"]
    available_node_features = None
    if features is not None:
        raise NotImplementedError("Only default features are supported for now")

    return g, {"edge": available_edge_features, "node": available_node_features}


@typechecked
def add_features(g: DbGraphType, features: Sequence[str] | None) -> Tuple[DbGraphType, FeatureDict]:
    if isinstance(g, nx.MultiDiGraph):
        if features:
            assert all(
                feat in supported_multidigraph_features() for feat in features
            ), "Unsupported feature attempted to be added"
        return add_multidigraph_features(g, features)
    else:
        raise ValueError(f"Unknown graph type {type(g)}")
