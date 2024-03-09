"""We need to load the graph from the file and construct the graph It is necessary to first load the edges in
dictionary (edge_id, edge_sequence) We need to load the nodes as overlap of the edges The nodes are defined as follows:

incoming edge -> outgoing edge and with the +/- defined as if incoming edge + use last kmer of the edge else from RC of
the edge if outgoing edge + use first kmer of the edge else from RC of the edge
"""
from pathlib import Path

import networkx as nx
from typeguard import typechecked

from graph.gfa_parser import SegmentDict, parse_gfa
from graph.rolling_hash import RollingHash

DbGraphType = nx.DiGraph | nx.MultiDiGraph


@typechecked
def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def verify_edge_overlaps(segments: SegmentDict, links: dict, k: int):
    for _, (inc_id, inc_sgn, out_id, out_sgn) in links.items():
        if inc_sgn == "+":
            inc = segments[inc_id][-k:]
            if out_sgn == "+":
                out = segments[out_id][:k]
                assert inc == out
            else:
                out = reverse_complement(segments[out_id][-k:])
                assert inc == out
        else:
            inc = reverse_complement(segments[inc_id][:k])
            if out_sgn == "+":
                out = segments[out_id][:k]
                assert inc == out
            else:
                out = reverse_complement(segments[out_id][-k:])
                assert inc == out


@typechecked
def construct_nx_multigraph(segments: SegmentDict, k: int) -> tuple[nx.MultiDiGraph, dict]:
    labeler = RollingHash(k=k)

    label_rc = {}
    g = nx.MultiDiGraph()
    for sid, attrs in segments.items():
        seq = attrs.pop("seq", None)
        if seq is None or seq == "*":
            raise ValueError(f"Invalid DNA sequence {seq}")

        kmer_start = seq[:k]
        kmer_end = seq[-k:]
        g.add_edge(kmer_start, kmer_end, key=sid, **attrs)

        kmer_rc_start = reverse_complement(kmer_end)
        kmer_rc_end = reverse_complement(kmer_start)

        sid_rc = labeler.hash(reverse_complement(seq))
        label_rc[sid] = sid_rc
        g.add_edge(kmer_rc_start, kmer_rc_end, key=sid_rc, **attrs)

    return g, label_rc


@typechecked
def construct_graph(graph_path: Path, k: int) -> tuple[DbGraphType, dict]:
    segments, _ = parse_gfa(path=graph_path, k=k)
    return construct_nx_multigraph(segments, k=k)
