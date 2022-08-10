from cgi import test
from graph.gfa_parser import parse_gfa


def test_parse_gfa(test_gfa_root):
    example_gfa = test_gfa_root / "example1.gfa"

    segments, links = parse_gfa(example_gfa, k=10)

    assert len(links) == 24
    assert len(segments) == 20