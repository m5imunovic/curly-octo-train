import filecmp
from pathlib import Path
from unittest import mock

from omegaconf import OmegaConf

import reference.reference_utils as ru
from reference.genome_generator import run as run_reference_step


def genome_generator_cfg(tmp_path):
    # TODO: replace with fixture loading the config and overwrite the path
    return OmegaConf.create(
        {
            "reference": {
                "species": {
                    "name": "test_species",
                    "chromosomes": {"chr1": 2000},
                    "gc_content": None,
                    "seed": 1,
                }
            },
            "paths": {"ref_dir": str(tmp_path / "curly_octo_train")},
        }
    )


@mock.patch("reference.reference_utils.save_chr_to_fasta")
def test_run_random_reference_saves_data(mock_save_chr_to_fasta, tmp_path):
    cfg = genome_generator_cfg(tmp_path=tmp_path)
    run_reference_step(cfg)

    expected_save_dir = ru.ref_chromosomes_path(ru.get_species_root(Path(cfg.paths.ref_dir), cfg.reference.species))
    assert mock_save_chr_to_fasta.call_count == 1
    assert expected_save_dir in mock_save_chr_to_fasta.call_args[0]
    assert all(chrX in mock_save_chr_to_fasta.call_args[0] for chrX in cfg.reference.species.chromosomes)

    # TODO: add test for real reference


def test_run_random_reference(tmp_path, test_species_cfg, test_data_reference):
    cfg = OmegaConf.create({"paths": {"ref_dir": str(tmp_path)}, "reference": {"species": test_species_cfg}})
    run_reference_step(cfg)
    save_dir = ru.ref_chromosomes_path(ru.get_species_root(Path(cfg.paths.ref_dir), cfg.reference.species))
    expected_dir = ru.ref_chromosomes_path(ru.get_species_root(test_data_reference, cfg.reference.species))

    assert filecmp.cmp(save_dir / "chr1.fasta", expected_dir / "chr1.fasta"), "The chr1.fasta file changed!"
    assert filecmp.cmp(save_dir / "chr2.fasta", expected_dir / "chr2.fasta"), "The chr2.fasta file changed!"
