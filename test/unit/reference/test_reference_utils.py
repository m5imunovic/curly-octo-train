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
