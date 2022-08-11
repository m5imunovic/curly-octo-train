from pathlib import Path
from typing import Dict

from icontract import require
from typeguard import typechecked


@require(lambda mult_info_path: mult_info_path.exists(), error=FileNotFoundError, description="Path does not exist.")
@typechecked
def parse_mult_info(mult_info_path: Path) -> Dict[str, int]:
    mult_info = {}
    with open(mult_info_path, 'r') as handle:
        for line in handle:
            name, multiplicity = line.strip().split()
            mult_info[name] = int(multiplicity)

    return mult_info
