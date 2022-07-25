import pytest
import tempfile
from collections import Counter
from pathlib import Path

from reference.genome_generator import get_random_chromosome, get_random_genome, save_genome_to_fasta


def test_get_random_chromosome():
    assert get_random_chromosome(10, seed=1) == 'GATCATGGTC'
    assert len(get_random_chromosome(10)) == 10


def test_get_random_genome_with_gc_content():
    gc_content_expected = 0.4
    genome = get_random_genome({'chr1': 1000000}, gc_content=gc_content_expected, seed=1)
    assert isinstance(genome, dict)

    chromosome_1 = genome['chr1']
    counter = Counter(chromosome_1)
    gc_content = (counter['G'] + counter['C']) / len(chromosome_1)
    assert pytest.approx(gc_content, 0.1) == gc_content_expected


def test_save_genome_to_fasta():
    with tempfile.TemporaryDirectory() as output_dir:
        path = Path(output_dir)
        save_genome_to_fasta(path, {'chr1': 'ACTGCTGATC'}, 'Test Genome')
        assert (path / 'chr1.fasta').exists()