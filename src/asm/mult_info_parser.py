from pathlib import Path
from typing import Dict, Set, Tuple

from icontract import require
from typeguard import typechecked


@require(lambda mult_info_path: mult_info_path.exists(), error=FileNotFoundError, description="Path does not exist.")
@typechecked
def parse_mult_info(mult_info_path: Path) -> Dict[str, int]:
    mult_info = {}
    with open(mult_info_path) as handle:
        for line in handle:
            name, multiplicity = line.strip().split()
            mult_info[name] = int(multiplicity)

    return mult_info


@typechecked
def partition_mult_info_edges(mult_info: Dict[str, int]) -> Tuple[Set, Set]:
    incorrect_edges = set()
    correct_edges = {key for key, value in mult_info.items() if value > 0 or incorrect_edges.add(key)}

    return correct_edges, incorrect_edges
