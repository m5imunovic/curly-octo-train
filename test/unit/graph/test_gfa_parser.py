import pytest

from graph.gfa_parser import parse_gfa


def test_parse_gfa(test_gfa_root):
    example_gfa = test_gfa_root / "example1.gfa"

    segments, links = parse_gfa(example_gfa, k=10)

    assert len(links) == 24
    assert len(segments) == 20

    segment1 = segments['1']
    assert 'seq' in segment1
    assert 'kc' in segment1
    assert 'ln' in segment1

    assert segment1['seq'] == '*'
    assert segment1['ln'] == 1
    assert pytest.approx(-0.22222222, rel=1e-6) == segment1['kc']