import logging
import re
from pathlib import Path

import networkx as nx
from typeguard import typechecked

logger = logging.getLogger(__name__)


@typechecked
def parse_dot(path: Path, k: int = 501) -> nx.MultiDiGraph:
    """Parse DBG from LJA graph.dot file"
    NOTE: this is thousand times slower than the custom implementation below.
    We don't need extra stuff as we know two possible line formats that appear in file

    Args:
        path (Path): Path to Dot file containing LaJolla DBG.
        k (int, optional): K-mer size used in JumboDBG stage. Defaults to 501.
    """
    g = nx.drawing.nx_pydot.read_dot(path)
    logger.info(f"Loaded graph with {g.number_of_nodes()} nodes and {g.number_of_edges} edges.")

    return g


@typechecked
def custom_parse_dot(path: Path, k: int = 501, verify=True) -> nx.MultiDiGraph:
    """Parses entries for edges
        "-2255" -> "-3481" [label="-2255.2 C 498(1 498)" color="black"]
        "2867" -> "665" [label="2867.1 T 66(16 1056)" color="black"]
    and nodes:
        "2867" [style=filled fillcolor="white"]
        "-665" [style=filled fillcolor="white"]

    Args:
        path (Path): Path to dot file
        k (int, optional): K-Mer size. Defaults to 501.
        verify: Checks that the nodes found by parsing node lines match nodes found in edge lines

    Raises:
        ValueError: If the file cannot be parsed

    Returns:
        nx.MultiDiGraph: De Bruijn assembly graph
    """
    nodes = set()
    edge_induced_nodes = set()
    g = nx.MultiDiGraph()
    with open(path) as f:
        for line in f:
            if line.startswith("digraph") or line.startswith("nodesep") or line.startswith("}"):
                continue
            l_id, l_meta = line.removesuffix("]\n").split(" [")
            if "->" in l_id:
                # edge entry
                start, stop = l_id.split(" -> ")
                start = start.strip('"')
                stop = stop.strip('"')
                # label and attributes
                l_label = l_meta.removesuffix(' color="black"').removeprefix("label=").strip('"')
                edge_id, _, trunc_size_cov = l_label.split(" ")
                trunc_size, cov = trunc_size_cov.split("(")
                cov = float(cov.strip(")"))
                edge_len = int(trunc_size) + k
                attrs = {"kc": cov, "ln": edge_len}
                # add edge to multidigraph
                g.add_edge(start, stop, key=edge_id, **attrs)
                if verify:
                    edge_induced_nodes.add(start)
                    edge_induced_nodes.add(stop)
            else:
                # node entry
                if verify:
                    nodes.add(l_id.strip('"'))
    if verify:
        assert nodes == edge_induced_nodes, "Nodes induced from edges don't match node list"

    return g


def custom_parse_dot2(path: Path, k: int = 501) -> nx.MultiDiGraph:
    """A bit more readable version of custom_parse_dot but also slower. We are trying to
    identify the rows with following format:
        "-2255" -> "-3481" [label="-2255.2 C 498(1 498)" color="black"]
        "2867" -> "665" [label="2867.1 T 66(16 1056)" color="black"]
    Args:
        path (Path): Path to dot file
        k (int, optional): K-Mer size. Defaults to 501.

    Raises:
        ValueError: If the file cannot be parsed

    Returns:
        nx.MultiDiGraph: De Bruijn assembly graph
    """
    g = nx.MultiDiGraph()
    with open(path) as f:
        for line in f:
            if line.startswith("digraph") or line.startswith("nodesep") or line.startswith("}"):
                continue
            if "->" not in line:
                continue

            # TODO: match this in two stages to make it more readable (edge part + attrs part)
            p = r"\"(?P<start>-?\d+)\" -> \"(?P<stop>-?\d+)\" \[label=\"(?P<label>[^\s]+) . (?P<trunc>\d+)\((?P<cov>[^\)]+).*"
            m = re.match(p, line)
            if m:
                entries = m.groupdict()
                start = entries["start"]
                stop = entries["stop"]
                label = entries["label"]
                ln = k + int(entries["trunc"])
                cov = float(entries["cov"])
                g.add_edge(start, stop, key=label, **{"cov": cov, "ln": ln})
            else:
                raise ValueError(f"Could not parse {line=}")

    return g
