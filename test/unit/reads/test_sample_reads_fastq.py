import filecmp
from pathlib import Path

import pytest

from reads.sample_reads_fastq import sample_reads_with_probability


@pytest.fixture(scope="module")
def test_data_profile_root(test_data_root) -> Path:
    return test_data_root / "pf"


def test_sample_reads_with_probability(tmpdir, test_data_profile_root):
    sampled_file = Path(tmpdir) / "sampled.fastq"
    reads_file = test_data_profile_root / "sample.fastq"

    selected = sample_reads_with_probability(reads_file, sampled_file, probability=0.25, seed=19)

    assert selected == 6, "No reads were selected!"
    expected_sampled_file = test_data_profile_root / "expected_sampled.fastq"
    assert filecmp.cmp(sampled_file, expected_sampled_file), "Sampled file does not match expected file!"
