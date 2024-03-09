from datetime import datetime
from graph.dot_parser import custom_parse_dot, custom_parse_dot2


def test_graph_dot_correctly(test_graph_root):
    graph_path = test_graph_root / "graph.dot"

    begin2 = datetime.now()
    g2 = custom_parse_dot(graph_path)
    duration = datetime.now() - begin2
    print(f"Parsing using custom took {duration}")

    begin3 = datetime.now()
    g3 = custom_parse_dot2(graph_path)
    duration = datetime.now() - begin3
    print(f"Parsing using custom2 took {duration}")

    assert g2.number_of_edges() == 9788
    assert g2.number_of_nodes() == 6722
    assert g3.number_of_edges() == 9788
    assert g3.number_of_nodes() == 6722