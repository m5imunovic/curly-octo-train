from reference.bed.bed_parser import parse_bed_file, merge_bed_regions


def test_parse_bed_file(test_data_reference):
    bed_file = test_data_reference / "bed" / "CHM13_VDJ.bed"
    bed_entries = parse_bed_file(bed_file)

    expected_keys = {"chr2", "chr7", "chr14", "chr22"}

    print(bed_entries)
    assert expected_keys == set(bed_entries)
    assert bed_entries["chr2"][0] == (88866370, 90790947)
    assert len(bed_entries["chr14"]) == 2


def test_merge_bed_regions(test_data_reference):
    bed_file = test_data_reference / "bed" / "hg002_chr10.bed"
    bed_entries = parse_bed_file(bed_file)
    merged_regions = merge_bed_regions(bed_entries)

    expected_keys = {"chr10_MATERNAL", "chr10_PATERNAL"}
    assert expected_keys == set(bed_entries)

    assert merged_regions["chr10_MATERNAL"] == (38477166, 44010545)
    assert merged_regions["chr10_PATERNAL"] == (38478676, 44593856)
