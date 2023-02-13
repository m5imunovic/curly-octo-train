"""Calculates the map of reverse complement IDs for reads in GFA files from LJA assembler.

Multiprocessing implementations take round third of the time of the single-threaded version.
"""
import multiprocessing as mp
from itertools import repeat
from pathlib import Path
from typing import Dict

from more_itertools import batched
from typeguard import typechecked

from graph.gfa_parser import parse_gfa
from graph.rolling_hash import RollingHash


def reverse_complement(seq: str) -> str:
    # TODO: unify this and the one in construct_graph.py
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


@typechecked
def get_rc_map(gfa_path: Path, k: int = 501) -> Dict:
    """Returns a dictionary mapping read IDs to their reverse complement."""
    # NOTE: this is awfully slow use one of the multiprocessing implementations below
    segments, _ = parse_gfa(gfa_path, k=k, skip_links=True)
    labels_rc = {}
    labeler = RollingHash(k=k)

    for sid, attrs in segments.items():
        seq = attrs.pop("seq", None)
        if seq is None or seq == "*":
            raise ValueError(f"Invalid DNA sequence {seq}")
        label = labeler.hash(reverse_complement(seq))
        labels_rc[sid] = label

    return labels_rc


def worker_batched(args, labeler):
    """The args argument is a chunk of input segments from GFA file Returns a dictionary mapping read IDs to their
    reverse complement."""
    labels_rc = {}
    for sid, attrs in args:
        seq = attrs.pop("seq", None)
        if seq is None or seq == "*":
            raise ValueError(f"Invalid DNA sequence {seq}")
        label = labeler.hash(reverse_complement(seq))
        labels_rc[sid] = label

    return labels_rc


def worker(args, labeler):
    """The args argument is a single item from input segments from GFA file Returns a dictionary mapping read ID to its
    reverse complement."""
    labels_rc = {}
    sid, attrs = args
    seq = attrs.pop("seq", None)
    if seq is None or seq == "*":
        raise ValueError(f"Invalid DNA sequence {seq}")
    label = labeler.hash(reverse_complement(seq))
    labels_rc[sid] = label

    return labels_rc


def get_rc_map_mp_pool_batch(gfa_path: Path, k: int = 501, threads: int = 8) -> Dict:
    """Returns a dictionary mapping read IDs to their reverse complement."""
    segments, _ = parse_gfa(gfa_path, k=k, skip_links=True)

    labeler = RollingHash(k=k)

    batch_size = len(segments) // threads
    print(f"Chunk size: {batch_size}")
    with mp.Pool(threads) as pool:
        processed_results = pool.starmap(worker_batched, zip(batched(segments.items(), batch_size), repeat(labeler)))

    labels_rc = {}
    for results in processed_results:
        labels_rc.update(results)

    return labels_rc


@typechecked
def get_rc_map_mp_pool(gfa_path: Path, k: int = 501, threads: int = 8) -> Dict:
    """Returns a dictionary mapping read IDs to their reverse complement."""
    segments, _ = parse_gfa(gfa_path, k=k, skip_links=True)

    labeler = RollingHash(k=k)

    chunk_size = threads
    print(f"Chunk size: {chunk_size}")
    with mp.Pool(threads) as pool:
        processed_results = pool.starmap(worker, zip(segments.items(), repeat(labeler)), chunksize=chunk_size)

    labels_rc = {}
    for results in processed_results:
        labels_rc.update(results)

    return labels_rc
