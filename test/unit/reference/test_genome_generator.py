import pytest
import tempfile
from collections import Counter
from pathlib import Path
from unittest import mock

from omegaconf import OmegaConf

from reference.genome_generator import run as run_reference_step
from reference.genome_generator import get_random_chromosome, get_random_genome, save_genome_to_fasta


def genome_generator_cfg(overwrite):
    return OmegaConf.create({
        'name': 'test_species',
        'overwrite': overwrite,
        'chromosomes': {
            'chr1': 2000
        },
        'gc_content': None,
        'seed': None
    })


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


@mock.patch('reference.genome_generator.save_chr_to_fasta')
@mock.patch('utils.path_helpers.get_ref_path', return_value=Path('/tmp/curly_octo_train'))
def test_overwrites_data(mock_get_ref_path, mock_save_chr_to_fasta):
    cfg = genome_generator_cfg(overwrite=True)
    run_reference_step(cfg)

    expected_save_dir = mock_get_ref_path.return_value / cfg.name / 'chromosomes'
    assert mock_save_chr_to_fasta.call_count == 1
    assert expected_save_dir in mock_save_chr_to_fasta.call_args[0]
    assert all(chrX in mock_save_chr_to_fasta.call_args[0] for chrX in cfg.chromosomes)
    run_reference_step(cfg)
    assert mock_save_chr_to_fasta.call_count == 2
    assert expected_save_dir in mock_save_chr_to_fasta.call_args[0]
    assert all(chrX in mock_save_chr_to_fasta.call_args[0] for chrX in cfg.chromosomes)

    cfg = genome_generator_cfg(overwrite=False)
    run_reference_step(cfg)
    assert mock_save_chr_to_fasta.call_count == 2

