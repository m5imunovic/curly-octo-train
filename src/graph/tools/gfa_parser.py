import logging
from pathlib import Path
from typing import Any

from typeguard import typechecked

from graph.tools.rolling_hash import RollingHash

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
        rh = RollingHash(k=k)
        for line in f:
            if line.startswith("S"):
                _, sid, seq, kc = line.strip().split()
                # LJA stores an id of forward and revese strand together separated by underline
                fw, rc = sid.split("_")
                kc = int(kc[len("KC:i:") :])
                ln = len(seq)

                hash_fw = rh.hash(seq)
                hash_rc = rh.hash(reverse_complement(seq))
                segments[fw] = {"hash": hash_fw, "kc": float(kc / (ln - k)), "ln": ln}
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
