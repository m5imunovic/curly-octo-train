from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from omegaconf import DictConfig, OmegaConf
from typeguard import typechecked


@dataclass
class ReferenceFiles:
    SPECIES_INFO: str = "species_info.json"


@typechecked
def is_random_reference(species: DictConfig) -> bool:
    return "seed" in species


@typechecked
def get_species_root(reference_root: Path, species: DictConfig) -> Path:
    if is_random_reference(species):
        return reference_root / species.name / f"S_{species.seed:04d}"
    else:
        return reference_root / species.name / species.release


@typechecked
def ref_chromosomes_path(species_root: Path) -> Path:
    return species_root / "chromosomes"


@typechecked
def check_reference_exists(references_root: Path, species_cfg: DictConfig) -> bool:
    # e.g. home_sapiens/chm13_v2 or random_species/S_0001
    # TODO: expand check to stat the files in species info
    species_root = get_species_root(references_root, species_cfg)
    return (species_root / ReferenceFiles.SPECIES_INFO).exists()


@typechecked
def create_species_info(species_cfg: DictConfig) -> dict:
    species_info = OmegaConf.to_object(species_cfg)
    species_info.update({"created": datetime.now().strftime("%Y-%m-%d %H:%M:%S")})
    return species_info


@typechecked
def save_chr_to_fasta(
    output_path: Path, chr_name: str, chr_seq: str, description: str = "", multiline: bool = True
) -> Path:
    file_path = output_path / f"{chr_name}.fasta"
    print(f"Saving {file_path}...")
    new_fasta = [SeqIO.SeqRecord(seq=Seq(chr_seq), id=chr_name, description=description)]
    if multiline:
        SeqIO.write(new_fasta, file_path, "fasta")
    else:
        # output entire DNA sequence in a single line
        new_fasta = [SeqIO.FastaIO.as_fasta_2line(record) for record in new_fasta]
        with open(file_path, "w") as handle:
            handle.writelines(new_fasta)

    return file_path


@typechecked
def save_genome_to_fasta(output_path: Path, genome: dict[str, str], description: str = "", multiline: bool = True):
    if not output_path.exists():
        output_path.mkdir(parents=True)
    for chr_name, chr_seq in genome.items():
        save_chr_to_fasta(output_path, chr_name, chr_seq, description, multiline)
