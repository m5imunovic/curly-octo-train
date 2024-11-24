import pytest

from graph.tools.gfa_parser import parse_gfa


def test_parse_lja_gfa(test_graph_root):
    example_gfa = test_graph_root / "gfa" / "final_dbg.gfa"

    segments, links = parse_gfa(example_gfa, k=501)

    assert len(links) == 0
    assert len(segments) == 10
    print(segments)

    assert "-1535.1" in segments.keys()
    assert "1871.2" in segments.keys()

    segment = segments["-1535.1"]
    assert "hash" in segment
    assert "kc" in segment
    assert "ln" in segment

    assert segment["ln"] == 1957
    assert pytest.approx(1456.0, rel=1e-6) == segment["kc"]
    assert segment["hash"] == "36757a92ce08eb96cf23bcf056ea239754309b75"
