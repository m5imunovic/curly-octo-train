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

    g = nx.MultiDiGraph()
    for sid, attrs in segments.items():
        assert attrs['seq'] is not '*'
        seq = attrs.pop('seq', None)
        kmer_start = seq[:k]
        kmer_end = seq[-k:]
        attrs.update({'sid': sid})
        g.add_edge(kmer_start, kmer_end, **attrs)

        kmer_rc_start = reverse_complement(kmer_end)
        kmer_rc_end = reverse_complement(kmer_start)

        attrs.update({'sid': '_' + sid})
        g.add_edge(kmer_rc_start, kmer_rc_end, **attrs)

    return g


@typechecked
def construct_nx_digraph(segments: SegmentDict, links: Dict[int, Tuple]) -> nx.DiGraph:
    g = nx.DiGraph()
    for sid, attrs in segments.items():
        attrs.pop('seq', None)
        g.add_node(sid, **attrs)
        g.add_node('_' + sid, **attrs)

    for (inc_id, inc_sgn, out_id, out_sgn) in links.values():
        vertex_from = inc_id if inc_sgn == '+' else '_' + inc_id
        vertex_to = out_id if out_sgn == '+' else '_' + out_id
        g.add_edge(vertex_from, vertex_to)

    return g


@typechecked
def construct_graph(cfg: DictConfig) -> Union[nx.DiGraph, nx.MultiDiGraph]:
    segments, links = parse_gfa(path=cfg.gfa_path, k=cfg.k)
    if cfg.graph_type == 'multigraph':
        return construct_nx_multigraph(segments, k=cfg.k)
    if cfg.graph_type == 'digraph':
        return construct_nx_digraph(segments, links)

    raise ValueError(f"Unknown graph type {cfg.graph_type}")
