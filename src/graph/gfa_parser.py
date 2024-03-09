import logging
from pathlib import Path
from typing import Any, Dict, Tuple

from typeguard import typechecked

SegmentDict = Dict[str, Dict[str, Any]]


logger = logging.getLogger(__name__)


@typechecked
def parse_gfa(path: Path, k: int = 501, skip_links=False) -> Tuple:
    """_summary_

    Args:
        path (Path):
            Path to GFA file containing LaJolla de Bruijn graph.
        k (int, optional):
            K-mer size used in JumboDBG stage. Defaults to 501.
        skip_links (bool, optional):
            Do not load links. Defaults to False.

    Returns:
        Tuple[SegmentDict, Dict]:
            Tuple of segments and links dictionaries.
    """
    segments = {}
    links = {}

    node_id = 0

    expected_cigar = str(k) + "M"

    with open(path) as f:
        version = f.readline().strip()
        assert "VN:Z:1.0" in version
        for line in f:
            if line.startswith("S"):
                _, sid, seq, kc = line.strip().split()
                kc = int(kc[len("KC:i:") :])
                ln = len(seq)
                segments[sid] = {"seq": seq, "kc": float(kc / (ln - k)), "ln": ln}
            if not skip_links:
                if line.startswith("L"):
                    _, inc_id, inc_sgn, out_id, out_sgn, cigar = line.strip().split()
                    assert cigar == expected_cigar
                    links[node_id] = (inc_id, inc_sgn, out_id, out_sgn)
                    node_id += 2

    logging.info(f"Loaded {len(segments)} segments")
    if not skip_links:
        logging.info(f"Loaded {len(links)} links")
    else:
        logging.info("Skipped links loading")

    return segments, links
