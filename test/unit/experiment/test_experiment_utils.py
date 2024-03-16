import re

import pytest

import experiment.experiment_utils as eu
from experiment.scenario_schema import Chromosome


def test_ensure_species_def_exists(test_cfg_root):
    # Assert that no exception is raised
    eu.ensure_species_def_exists(test_cfg_root, ["test_species.yaml"])


def test_ensure_species_def_exists_missing_file(test_cfg_root):
    with pytest.raises(
        AssertionError,
        match=r"Species definition does not exist at .*/reference/species/non_existing_species_def.yaml",
    ):
        eu.ensure_species_def_exists(test_cfg_root, ["non_existing_species_def.yaml"])


def test_ensure_species_def_exists_invalid_extension(test_cfg_root):
    # Call the function with the species definition file with invalid extension
    with pytest.raises(AssertionError, match="Species definition file must have a .yaml extension: species_def.txt"):
        eu.ensure_species_def_exists(test_cfg_root, ["species_def.txt"])


def test_get_real_read_path(test_species_reads_root):
    reads = eu.get_read_paths(test_species_reads_root, chromosomes=[Chromosome(name="chr1")])
    expected_file1 = test_species_reads_root / "chr1" / "chr1_part1.fastq"
    expected_file2 = test_species_reads_root / "chr1" / "chr1_part2.fastq"
    assert len(reads) == 2
    assert expected_file1 in reads
    assert expected_file2 in reads

    reads = eu.get_read_paths(test_species_reads_root, chromosomes=[Chromosome(name="chr3")])
    assert not reads, "Unexpected reads found for chromosome chr3"


def test_get_sequencing_seeds():
    subset, count, init_seed = "train", 5, 1
    expected_seeds = [10001, 10002, 10003, 10004, 10005]
    assert eu.get_sequencing_seeds(subset, count, init_seed) == expected_seeds

    subset, count, init_seed = "val", 3, 8765
    expected_seeds = [28765, 28766, 28767]
    assert eu.get_sequencing_seeds(subset, count, init_seed) == expected_seeds

    subset, count, init_seed = "test", 1, 4321
    expected_seeds = [34321]
    assert eu.get_sequencing_seeds(subset, count, init_seed) == expected_seeds


def test_get_sequencing_seeds_erroneous_input():
    with pytest.raises(
        AssertionError, match=re.escape("Unknown subset: unknown, expected one of ['train', 'val', 'test']")
    ):
        eu.get_sequencing_seeds(subset="unknown", count=1, init_seed=1)
    with pytest.raises(AssertionError, match="Sum of initial seed and count must be less than 10000"):
        eu.get_sequencing_seeds(subset="test", count=100, init_seed=9900)
    with pytest.raises(AssertionError, match="Initial seed must be between 1 and 9999"):
        eu.get_sequencing_seeds(subset="test", count=100, init_seed=10000)


def test_get_sequencing_probability():
    none_probability = None
    probabilities = eu.get_sequencing_probabilities(none_probability) 
    assert probabilities == [1.0]

    probability_range = "0.92:0.92:0.01"
    probabilities = eu.get_sequencing_probabilities(probability_range)
    assert probabilities == pytest.approx([0.92])

    probability_range = "0.92:0.98:0.02"
    probabilities = eu.get_sequencing_probabilities(probability_range)
    assert probabilities == pytest.approx([0.92, 0.94, 0.96, 0.98])

    probability_range = "0.98:1.00:0.02"
    probabilities = eu.get_sequencing_probabilities(probability_range)
    assert probabilities == pytest.approx([0.98, 1.00])
    
    probability_range = "0.98:1.00:0.25"
    probabilities = eu.get_sequencing_probabilities(probability_range)
    assert probabilities == pytest.approx([0.98])

