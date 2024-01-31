from collections import Counter

import numpy as np
import pytest

from reference.random_genome import get_random_chromosome, get_random_reference


def test_get_random_chromosome():
    rng = np.random.default_rng(42)
    sequence_length = 100
    chromosome = get_random_chromosome(sequence_length, rng)

    assert len(chromosome) == sequence_length
    assert set(chromosome) == set("ACGT")


def test_get_random_genome_with_gc_content():
    gc_content_expected = 0.4
    genome = get_random_reference({"chr1": 1000000}, seed=1, gc_content=gc_content_expected)
    assert isinstance(genome, dict)

    chromosome_1 = genome["chr1"]
    counter = Counter(chromosome_1)
    gc_content = (counter["G"] + counter["C"]) / len(chromosome_1)
    assert pytest.approx(gc_content, 0.1) == gc_content_expected


def test_get_random_genome():
    chromosomes = {"chr1": 100, "chr2": 100, "chr3": 150}
    seed = 42
    genome = get_random_reference(chromosomes, seed)

    assert isinstance(genome, dict)
    assert len(genome) == len(chromosomes)

    for chromosome, sequence in genome.items():
        assert len(sequence) == chromosomes[chromosome]
        assert set(sequence) == set("ACGT")

    assert genome["chr1"] != genome["chr2"]
    assert genome["chr1"][:5] == "GATCA"
    assert genome["chr2"][:5] == "CCGAG"
    assert genome["chr3"][:5] == "CAAAG"
