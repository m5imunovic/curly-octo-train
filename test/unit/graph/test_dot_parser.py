from datetime import datetime

from graph.dot_parser import custom_parse_dot, custom_parse_dot2


def test_graph_dot_correctly(test_dot_root):
    graph_path = test_dot_root / "example1.dot"

    begin2 = datetime.now()
    g2 = custom_parse_dot(graph_path)
    duration = datetime.now() - begin2
    print(f"Parsing using custom took {duration}")

    begin3 = datetime.now()
    g3 = custom_parse_dot2(graph_path)
    duration = datetime.now() - begin3
    print(f"Parsing using custom2 took {duration}")

    assert g2.number_of_nodes() == 3502
    assert g2.number_of_edges() == 5116
    assert g3.number_of_nodes() == 3502
    assert g3.number_of_edges() == 5116
