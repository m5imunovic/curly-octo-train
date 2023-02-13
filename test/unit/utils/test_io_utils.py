import tempfile
from pathlib import Path

from utils.io_utils import compose_cmd_params, get_read_files


def test_compose_cmd_params():
    params_empty = {}
    assert compose_cmd_params(params_empty) == ""

    params_short_no_val = {"short": {"a": ""}}
    assert compose_cmd_params(params_short_no_val) == "-a"

    params_long_no_val = {"long": {"a": ""}}
    assert compose_cmd_params(params_long_no_val) == "--a"

    params_short = {"short": {"a": "1", "b": "2"}}
    assert compose_cmd_params(params_short) == "-a 1 -b 2"

    params_long = {"long": {"a": "1", "b": "2"}}
    assert compose_cmd_params(params_long) == "--a 1 --b 2"

    params = {"short": {"a": "1", "b": "2"}, "long": {"c": "3", "d": "4"}, "append": "-e 5"}
    assert compose_cmd_params(params) == "-a 1 -b 2 --c 3 --d 4 -e 5"

    params = {
        "short": {"a": "1", "b": "2"},
        "long": {"with_underscore": "3", "with-dash": "4"},
        "append": "whatever the argument is",
    }

    assert compose_cmd_params(params) == "-a 1 -b 2 --with_underscore 3 --with-dash 4 whatever the argument is"


def test_get_read_files(test_reads_root):
    expected_read_files = ["dummy.fa", "dummy.fastq"]
    read_files = get_read_files(test_reads_root / "dummy", suffix=[".fasta", ".fa"])
    assert all([f.name in expected_read_files for f in read_files])

    expected_read_files = ["dummy.fastq"]
    read_files = get_read_files(test_reads_root / "dummy", suffix=[".fastq"], override=True)
    assert all([f.name in expected_read_files for f in read_files])

    with tempfile.NamedTemporaryFile(suffix=".fastq") as handle:
        expected_read_files = [Path(handle.name).name]
        read_files = get_read_files(Path(handle.name), override=False)
        assert len(read_files) == 1
        assert all([f.name in expected_read_files for f in read_files])
