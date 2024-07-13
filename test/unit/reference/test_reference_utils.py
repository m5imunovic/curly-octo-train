import tempfile
from pathlib import Path

from omegaconf import DictConfig

import reference.reference_utils as ru


def test_get_species_root():
    reference_root = Path("/reference")
    species = DictConfig({"name": "homo_sapiens", "release": "GRCh38"})
    assert ru.get_species_root(reference_root, species) == Path("/reference/homo_sapiens/GRCh38")

    species = DictConfig({"name": "random_species", "seed": 42})
    assert ru.get_species_root(reference_root, species) == Path("/reference/random_species/S_0042")


def test_reference_exists(test_data_reference, test_species_cfg, tmp_path):
    assert ru.check_reference_exists(test_data_reference, test_species_cfg)
    assert not ru.check_reference_exists(tmp_path, test_species_cfg)


def test_save_genome_to_fasta():
    with tempfile.TemporaryDirectory() as output_dir:
        path = Path(output_dir)
        ru.save_genome_to_fasta(path, {"chr1": "ACTGCTGATC"}, "Test Genome")
        assert (path / "chr1.fasta").exists()


def test_query_chr_paths(test_data_reference, test_species_cfg):
    species_root = ru.get_species_root(test_data_reference, test_species_cfg)
    chromosomes_path = ru.ref_chromosomes_path(species_root)
    chromosome_paths = ru.query_chr_paths(chromosomes_path)
    assert len(chromosome_paths) == 2
    assert "chr1" in chromosome_paths
    assert "chr2" in chromosome_paths
    assert chromosome_paths["chr1"].name == "chr1.fasta"


def test_extract_chromosome_range():
    with tempfile.TemporaryDirectory() as output_dir:
        path = Path(output_dir)
        ru.save_genome_to_fasta(path, {"chr1": "ACTGCTGATC"}, "Test Genome")
        chr1_path = path / "chr1.fasta"
        assert chr1_path.exists()
        chr_subref = ru.extract_chromosome_range(chr1_path, subrange=(2, 5))
        assert "chr1" in chr_subref
        assert chr_subref["chr1"] == "TGC"
