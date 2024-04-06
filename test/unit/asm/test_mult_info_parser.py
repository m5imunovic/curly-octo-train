import pytest

from asm.mult_info_parser import parse_mult_info, partition_mult_info_edges


def test_parse_mult_info(test_dot_root):
    mult_info_path = test_dot_root / "example1_mult.info"
    mult_info = parse_mult_info(mult_info_path)
    assert len(mult_info) == 5116
    assert "5.11" in mult_info
    assert "-9.10" in mult_info
    assert mult_info["5.11"] == 1
    assert mult_info["-9.10"] == 0


def test_parse_mult_info_not_exists(test_graph_root):
    mult_info_path = test_graph_root / "example_not_exists_mult.info"
    with pytest.raises(FileNotFoundError):
        parse_mult_info(mult_info_path)


def test_partition_mult_info_edges():
    mult_info = {"A": 1, "B": 0, "C": 5, "D": 0, "E": 10}

    correct_edges, incorrect_edges = partition_mult_info_edges(mult_info)
    assert correct_edges == {"A", "C", "E"}
    assert incorrect_edges == {"B", "D"}
