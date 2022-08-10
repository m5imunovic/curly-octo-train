from pathlib import Path
from typing import Dict, Tuple, Set


from graph.rolling_hash import RollingHash


def simple_hash(old_id: str, seq: str, k: int = 501):
    assert len(seq) > k
    new_id = old_id[:-2]
    new_id += '0' if old_id[-2] == '1' else '1'
    new_id += seq[k:k+1]

    return new_id



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


def parse_gfa2(path: Path, k=501) -> Tuple[Set, Set]:
    nodes = set()
    nodes_impl = set()
    edges = set()
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('S'):
                parsed = line.strip().split()
                sid = parsed[1]
                seq = parsed[2]
                kc = parsed[3][len("KC:i:")]
                nodes.add(sid)
                rh = RollingHash(k=k)
                rh_ = simple_hash(sid, seq, k=k)
                h, h_rc = rh.hash(seq, pos=0)
                if h == sid:
                    nodes_impl.add(h_rc)
                elif h_rc == sid:
                    nodes_impl.add(h)


            if line.startswith('L'):
                parsed = line.strip().split()
                from_id = parsed[1]
                to_id = parsed[3]
                edges.add((from_id, to_id))

    
    for from_id, to_id in edges:
        assert from_id in nodes
        assert to_id in nodes


    return nodes, edges