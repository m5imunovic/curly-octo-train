import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from asm import la_jolla
from asm.assembler import assembler_factory
from asm.mult_info_parser import parse_mult_info


@pytest.fixture(scope="module")
def fake_vendor_root():
    # TODO: use pyfakefs for this
    with tempfile.TemporaryDirectory() as tmp_dir:
        vendor_dir = Path(tmp_dir) / "vendor"
        vendor_dir.mkdir(exist_ok=True)
        assembler_root = vendor_dir / "LJA"
        assembler_root.mkdir(exist_ok=True)
        # LaJolla expects vendor_dir as input argument and that LJA folder exists in there
        yield vendor_dir


def test_la_jolla_construct_exe_cmd(fake_vendor_root, test_la_jolla_cfg, tmpdir):
    lja = la_jolla.LaJolla(test_la_jolla_cfg, fake_vendor_root)
    output_path = Path(tmpdir) / "fake/output"
    reads = [Path("/fake/reads/chr1/0_0001.fastq"), Path("/fake/reads/chr1/1_0001.fastq")]
    genome = [Path("/fake/reference/chromosomes/chr1.fasta")]
    with patch("json.dump") as mock_dump, patch("builtins.open") as mock_open:
        cmd, eval_cmd_path = lja._construct_exec_cmd(genome=genome, reads=reads, output_path=output_path)

    # call jumboDBG
    assert len(cmd) == 1
    assert str(fake_vendor_root / "LJA/bin/mlgraph") in cmd[0]
    assert "-k 501" in cmd[0]
    assert "--threads 15" in cmd[0]
    assert "--compress" in cmd[0]
    assert f"--reference {genome[0]}" in cmd[0]
    assert f"--reads {reads[0]}" in cmd[0]
    assert f"--reads {reads[1]}" in cmd[0]

    # prepares the commands for later evaluation
    # call lja
    # call align_and_print twice
    assert output_path / test_la_jolla_cfg["full_asm"] / "full_asm.json" == eval_cmd_path
    assert "full_asm.json" in str(mock_open.call_args_list[0])
    assert "lja" in mock_dump.call_args[0][0]
    assert "align_and_print" in mock_dump.call_args[0][0]


@pytest.mark.xfail(reason="Need to update the files")
def test_la_jolla_produces_expected_output(
    test_la_jolla_cfg, test_species_genome, test_data_reads_fq, test_data_assemblies, tmpdir
):
    assembler = assembler_factory("LJA", test_la_jolla_cfg)

    output_path = Path(tmpdir)
    # with patch("json.dump"), patch("builtins.open"):
    assembler(genome=test_species_genome, reads=test_data_reads_fq, output_path=output_path)

    assert output_path.exists()

    assert (output_path / "graph.dot").exists()
    assert (output_path / "graph.gfa").exists()
    assert (output_path / "mult.info").exists()
    assert (output_path / "full_asm" / "full_asm.json").exists()
    expected_mult_info = parse_mult_info(test_data_assemblies / "mult.info")
    mult_info = parse_mult_info(output_path / "mult.info")
    assert sorted(mult_info.values()) == sorted(expected_mult_info.values())
    # after switching to the number based IDs keys are not same across the runs
    assert len(mult_info.keys()) == len(expected_mult_info.keys())
