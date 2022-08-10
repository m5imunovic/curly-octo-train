from pathlib import Path
from typing import Dict, Tuple, Set
from typeguard import typechecked


@typechecked
def parse_gfa(path: Path, k=501) -> Tuple[Dict, Dict]:
    segments = {}
    links = {}

    node_id = 0

    expected_cigar = str(k) + 'M'

    with open(path, 'r') as f:
        version = f.readline().strip()
        assert 'VN:Z:1.0' in version
        for line in f:
            if line.startswith('S'):
                _, sid, seq, kc = line.strip().split()
                kc = kc[len("KC:i:")]
                segments[sid] = seq
            
            if line.startswith('L'):
                _, inc_id, inc_sgn, out_id, out_sgn, cigar = line.strip().split()
                assert cigar == expected_cigar
                links[node_id] = (inc_id, inc_sgn, out_id, out_sgn)
                node_id += 2

    print(f'Loaded {len(segments)} segments')
    print(f'Loaded {len(links)} links')

    return segments, links
