"""Adds additional features to the graph based on config."""
from typing import Optional, Sequence, Union

import networkx as nx
from typeguard import typechecked

DbGraphType = nx.MultiDiGraph
FeatureDict = Union[dict[str, Optional[list]], all]


def supported_multidigraph_features() -> set[str]:
    return {"ln", "kc"}


@typechecked
def add_multidigraph_features(
    g: nx.MultiDiGraph, features: Sequence[str] | None
) -> tuple[nx.MultiDiGraph, FeatureDict]:
    available_edge_features = ["kc", "ln"]
    available_node_features = None
    if features is not None:
        raise NotImplementedError("Only default features are supported for now")

    return g, {"edge": available_edge_features, "node": available_node_features}


@typechecked
def add_features(g: DbGraphType, features: Sequence[str] | None) -> tuple[DbGraphType, FeatureDict]:
    if isinstance(g, nx.MultiDiGraph):
        if features:
            assert all(
                feat in supported_multidigraph_features() for feat in features
            ), "Unsupported feature attempted to be added"
        return add_multidigraph_features(g, features)
    else:
        raise ValueError(f"Unknown graph type {type(g)}")


def add_mult_info_features(g: DbGraphType, mult_info: dict[str, int]) -> DbGraphType:
    # We use a small trick here. As we name this feature 'y' it automatically gets assigned to the y
    # attribute of the exported pytorch geometric graph and can then be used for learning supervision
    if isinstance(g, nx.MultiDiGraph):
        for edge in g.edges(data=True, keys=True):
            key = edge[2]
            g.edges[edge[0], edge[1], key]["y"] = mult_info[key]
    else:
        raise TypeError(f"Graph type {type(g)} not supported")
    return g
