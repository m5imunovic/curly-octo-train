import shutil

from Bio import SeqIO
from omegaconf import OmegaConf

import reference.reference_utils as ru
from reference.bed_genome_generator import run as generate_bed_genome


def test_generate_bed_reference(tmp_path, test_real_species_cfg, test_data_reference):
    bed_file_path = test_data_reference / "bed" / "test_species.bed"
    cfg = OmegaConf.create(
        {
            "paths": {"ref_dir": str(test_data_reference)},
            "reference": {"species": test_real_species_cfg, "bed": {"url": f"file://{bed_file_path}"}},
        }
    )

    generate_bed_genome(cfg, neighborhood=0, gap=100)
    expected_bed_species = OmegaConf.create({"name": "test_species", "release": "bed_S_0001"})
    bed_reference_dir = ru.get_species_root(test_data_reference, expected_bed_species)

    chr2_fasta_path = bed_reference_dir / "chromosomes" / "chr2.fasta"
    assert chr2_fasta_path.exists()
    entry = list(SeqIO.parse(chr2_fasta_path, "fasta"))[0]
    assert len(entry.seq) == 60
    assert entry.id == "chr2"
    shutil.rmtree(bed_reference_dir)
