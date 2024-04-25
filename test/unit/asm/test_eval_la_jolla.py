from asm.eval_la_jolla import parse_alignments_entry


def alignment_entry():
    entry = (
        ">S1_4\n"
        "-717.13 -1373.11 836.11 -1194.12 1913.13 -1475.11 -1459.10 92.21 -1419.11 1348.13 310.13 -223.10 -684.11 -1123.10 -634.11 245.13 -385.12 316.11 1344.12 -1898.11 -604.10 -1645.10 149.12 -1452.12 273.12 1729.13 -326.11 -1365.10 -922.13 74.13 386.12 -1174.20"
    )
    yield from entry.split("\n")


def test_parse_alignments_entry():
    entry = alignment_entry()
    for read_line in entry:
        read_id, edge_ids = parse_alignments_entry(read_line, next(entry))
        assert read_id == "S1_4"
        assert edge_ids[0] == "-717.13"
        assert edge_ids[2] == "836.11"
        assert edge_ids[-1] == "-1174.20"
