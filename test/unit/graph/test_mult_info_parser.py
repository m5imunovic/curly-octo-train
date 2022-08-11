import pytest

from graph.mult_info_parser import parse_mult_info


def test_parse_mult_info(test_graph_root):
    mult_info_path = test_graph_root / 'example_mult.info'
    mult_info = parse_mult_info(mult_info_path)
    assert len(mult_info) == 10
    assert '3366311032649037625441063611614947082921T' in mult_info


def test_parse_mult_info_not_exists(test_graph_root):
    mult_info_path = test_graph_root / 'example_not_exists_mult.info'
    with pytest.raises(FileNotFoundError):
        parse_mult_info(mult_info_path)