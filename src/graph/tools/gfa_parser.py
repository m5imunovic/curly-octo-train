import hashlib
import json
import logging
from pathlib import Path
from typing import Any

from typeguard import typechecked

SegmentDict = dict[str, dict[str, Any]]

logger = logging.getLogger(__name__)


@typechecked
def reverse_complement(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


@typechecked
def parse_gfa(path: Path, k: int = 501, skip_links: bool = True) -> tuple:
    """Custom parser for gfa produced by LJA assembler. Edgse are stored as segments and nodes are stored as overlap
    (length = 501) of two edges.

    Args:
        path (Path):
            Path to GFA file containing LaJolla de Bruijn graph.
        k (int, optional):
            K-mer size used in JumboDBG stage. Defaults to 501.
        skip_links (bool, optional):
            Do not load links. Defaults to True.

    Returns:
        tuple[SegmentDict, dict]:
            Tuple of segments and links dictionaries.
    """
    segments = {}
    links = {}

    str(k) + "M"

    with open(path) as f:
        version = f.readline().strip()
        assert "VN:Z:1.0" in version
        for line in f:
            if line.startswith("S"):
                _, sid, seq, kc = line.strip().split()
                # LJA stores an id of forward and revese strand together separated by underline
                kc = int(kc[len("KC:i:") :])
                ln = len(seq)

                split_id = sid.split("_")
                fw = split_id[0]
                hash_fw = hashlib.sha1(seq.encode("utf-8"), usedforsecurity=False).hexdigest()
                segments[fw] = {"hash": hash_fw, "kc": float(kc / (ln - k)), "ln": ln}

                has_rc = len(split_id) == 2
                if has_rc:
                    rc = split_id[1]
                    hash_rc = hashlib.sha1(seq.encode("utf-8"), usedforsecurity=False).hexdigest()
                    segments[rc] = {"hash": hash_rc, "kc": float(kc / (ln - k)), "ln": ln}
            if not skip_links:
                raise NotImplementedError("This functionality is not ported to new format yet")
    #                 if line.startswith("L"):
    #                     _, inc_id, inc_sgn, out_id, out_sgn, cigar = line.strip().split()
    #                     assert cigar == expected_cigar
    #                     links[node_id] = (inc_id, inc_sgn, out_id, out_sgn)
    #                     node_id += 2

    logger.info(f"Loaded {len(segments)} segments")
    if not skip_links:
        logger.info(f"Loaded {len(links)} links")
    else:
        logger.info("Skipped links loading")

    return segments, links


@typechecked
def save_gfa_hashmap(gfa: dict, output_path: Path):
    """Saves the result of `parser_gfa` into file for use in evaluation phase."""

    hashmap = {}
    for seg_id, metadata in gfa.items():
        hashmap[metadata["hash"]] = seg_id

    with open(output_path, "w") as f:
        json.dump(hashmap, f)
