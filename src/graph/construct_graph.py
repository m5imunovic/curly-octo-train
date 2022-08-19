"""
We need to load the graph from the file and construct the graph
It is necessary to first load the edges in dictionary (edge_id, edge_sequence)
We need to load the nodes as overlap of the edges
The nodes are defined as follows:
incoming edge -> outgoing edge
and with the +/- defined as
if incoming edge + use last kmer of the edge else from RC of the edge
if outgoing edge + use first kmer of the edge else from RC of the edge
"""
from typing import Dict, Tuple, Union

import networkx as nx
from omegaconf import DictConfig
from typeguard import typechecked

from graph.gfa_parser import parse_gfa, SegmentDict
from graph.rolling_hash import RollingHash

DbGraphType = Union[nx.DiGraph, nx.MultiDiGraph]


@typechecked
def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def verify_edge_overlaps(segments: SegmentDict, links: dict, k: int):
    for _, (inc_id, inc_sgn, out_id, out_sgn) in links.items():
        if inc_sgn == '+':
            inc = segments[inc_id][-k:]
            if out_sgn == '+':
                out = segments[out_id][:k]
                assert inc == out
            else:
                out = reverse_complement(segments[out_id][-k:])
                assert inc == out
        else:
            inc = reverse_complement(segments[inc_id][:k])
            if out_sgn == '+':
                out = segments[out_id][:k]
                assert inc == out
            else:
                out = reverse_complement(segments[out_id][-k:])
                assert inc == out


@typechecked
def construct_nx_multigraph(segments: SegmentDict, k: int) -> nx.MultiDiGraph:
    labeler = RollingHash(k=k)

    g = nx.MultiDiGraph()
    for sid, attrs in segments.items():
        seq = attrs.pop('seq', None)
        if seq is None or seq == '*':
            raise ValueError(f'Invalid DNA sequence {seq}')

        kmer_start = seq[:k]
        kmer_end = seq[-k:]
        attrs.update({'sid': sid})
        g.add_edge(kmer_start, kmer_end, **attrs)

        kmer_rc_start = reverse_complement(kmer_end)
        kmer_rc_end = reverse_complement(kmer_start)

        sid_rc = labeler.hash(reverse_complement(seq))
        attrs.update({'sid': sid_rc})
        g.add_edge(kmer_rc_start, kmer_rc_end, **attrs)

    return g


@typechecked
def construct_nx_digraph(segments: SegmentDict, links: Dict[int, Tuple], k: int) -> nx.DiGraph:
    g = nx.DiGraph()
    labels_rc = {}
    labeler = RollingHash(k=k)

    for sid, attrs in segments.items():
        seq = attrs.pop('seq', None)
        if seq is None or seq == '*':
            raise ValueError(f'Invalid DNA sequence {seq}')
        g.add_node(sid, **attrs)
        label = labeler.hash(reverse_complement(seq))
        labels_rc[sid] = label
        g.add_node(label, **attrs)

    for (inc_id, inc_sgn, out_id, out_sgn) in links.values():
        vertex_from = inc_id if inc_sgn == '+' else labels_rc[inc_id]
        vertex_to = out_id if out_sgn == '+' else labels_rc[out_id]
        g.add_edge(vertex_from, vertex_to)

    return g


@typechecked
def construct_graph(cfg: DictConfig) -> DbGraphType:
    segments, links = parse_gfa(path=cfg.gfa_path, k=cfg.k)
    if cfg.graph_type == 'multigraph':
        return construct_nx_multigraph(segments, k=cfg.k)
    if cfg.graph_type == 'digraph':
        return construct_nx_digraph(segments, links, k=cfg.k)

    raise ValueError(f"Unknown graph type {cfg.graph_type}")
