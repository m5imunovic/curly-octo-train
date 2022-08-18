from graph.construct_graph import reverse_complement
from graph.gfa_parser import parse_gfa
from graph.mult_info_parser import parse_mult_info
from graph.rolling_hash import RollingHash


def test_generate_rolling_hash(test_gfa_root):
    rh = RollingHash()
    example_gfa = test_gfa_root / 'lja_graph.gfa'
    segments, _ = parse_gfa(example_gfa, k=rh.k)

    for sid in segments:
        seq = segments[sid]['seq']
        hash = rh.hash(seq)

        assert sid == hash
    

def test_generate_rolling_hash_pair(test_gfa_root, test_graph_root):
    rh = RollingHash()

    example_gfa_path = test_gfa_root / 'lja_graph.gfa'
    segments, _ = parse_gfa(example_gfa_path, k=rh.k)

    hashes = set()
    hashes_ = set()
    for sid in segments:
        seq = segments[sid]['seq']
        hash = rh.hash(seq, pos=0)
        hashes.add(hash)
        hash_ = rh.hash(reverse_complement(seq), pos=0)
        hashes_.add(hash_)


    mult_info_path = test_graph_root / 'lja_mult.info'
    mult_info = parse_mult_info(mult_info_path)

    assert len(mult_info) == len(hashes) + len(hashes_)

    diff = set(mult_info.keys()) - hashes - hashes_
    assert len(diff) == 0