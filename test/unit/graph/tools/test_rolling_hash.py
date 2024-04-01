import json

from graph.tools.rolling_hash import RollingHash


def test_generate_rolling_hash(test_graph_root):
    rh = RollingHash()
    rolling = test_graph_root / "gfa" / "rolling.json"
    with open(rolling) as f:
        data = json.load(f)
        assert data["hash"] == rh.hash(data["seq"])
