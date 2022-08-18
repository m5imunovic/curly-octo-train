from cgi import test
import pytest

from graph.gfa_parser import parse_gfa


def test_parse_gfa(test_gfa_root):
    example_gfa = test_gfa_root / "example1.gfa"

    segments, links = parse_gfa(example_gfa, k=10)

    assert len(links) == 24
    assert len(segments) == 20

    segment = segments['1']
    assert 'seq' in segment
    assert 'kc' in segment
    assert 'ln' in segment

    assert segment['seq'] == '*'
    assert segment['ln'] == 1
    assert pytest.approx(-0.22222222, rel=1e-6) == segment['kc']


# TODO: parametrize test in order to merge these two into a single function
def test_parse_lja_gfa(test_gfa_root):
    example_gfa = test_gfa_root / 'lja_graph.gfa'

    segments, links = parse_gfa(example_gfa, k=501)

    assert len(links) == 66
    assert len(segments) == 80

    segment = segments['1227474217660589217047698276682431851491C']
    assert 'seq' in segment
    assert 'kc' in segment
    assert 'ln' in segment

    assert segment['ln'] == 1431
    assert segment['seq'].startswith('ATGCTGTCGAGCACTATGAC')
    assert pytest.approx(0.0096774193, rel=1e-6) == segment['kc']