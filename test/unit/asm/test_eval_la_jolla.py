from asm.eval_la_jolla import parse_alignments_entry


def alignment_entry():
    entry = (
        ">S1_4\n"
        "--597.2_762.1 597.1_-954.1 954.1_1618.2 -1471.2_-1618.1 -1117.1_-1471.1 -5.1_-1117.2 -5.2_-1300.2 1300.1_-1323.1 -825.1_1323.2 -825.2_1752.1 -802.1_-1752.1 -611.1_-802.1 --22.1_-611.1 22.2_990.1 -990.2_1134.2 --714.1_-1134.1 714.2_1225.1 -1225.1_1750.1 -434.1_-1750.1 --277.1_-434.1 277.1_-1888.1 --531.1_1888.1 531.1_-636.1 -536.1_636.1 -536.1_1515.1 -408.2_-1515.1 -408.1_752.1 -164.2_-752.1 -164.1_-1325.1 -665.1_1325.2 -649.1_-665.1 -649.1_-1307.1 1307.1_1917.1 "
    )
    yield from entry.split("\n")


def test_parse_alignments_entry():
    entry = alignment_entry()
    for read_line in entry:
        read_id, edge_ids_with_rc, is_rc = parse_alignments_entry(read_line, next(entry))
        assert read_id == "S1_4"
        assert edge_ids_with_rc[0] == "-597.2_762.1"
        assert edge_ids_with_rc[1] == "597.1_-954.1"
        assert edge_ids_with_rc[-1] == "1307.1_1917.1"
        assert is_rc[0] is True
        assert is_rc[1] is False
        assert is_rc[-1] is False
