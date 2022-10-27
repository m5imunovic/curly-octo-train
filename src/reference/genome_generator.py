from pathlib import Path
from typing import Dict, Optional, Sequence

import hydra
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from omegaconf import DictConfig
from typeguard import typechecked

from utils import path_helpers as ph


@typechecked
def get_random_chromosome(sequence_length: int, base_p: Optional[Sequence[float]] = None, seed: Optional[int] = None) -> str:
    if seed:
        np.random.seed(seed)
    rng = np.random.default_rng(12345)
    bases = rng.choice(['A', 'C', 'G', 'T'], size=sequence_length, p=base_p)
    return ''.join(bases)


@typechecked
def get_random_genome(chromosomes: Dict[str, int], gc_content: Optional[float] = None, seed: Optional[int] = None) -> Dict[str, str]:

    base_p = [(1-gc_content)/2, gc_content/2, gc_content/2, (1-gc_content)/2] if gc_content else None

    genome = {}
    for chr_name, chr_length in chromosomes.items():
        chr_seq = get_random_chromosome(chr_length, base_p, seed)
        genome[chr_name] = chr_seq

    return genome


@typechecked
def save_chr_to_fasta(output_path: Path, chr_name: str, chr_seq: str, description: str = '', multiline: bool = True) -> Path:
    file_path = output_path / f'{chr_name}.fasta'
    new_fasta = [SeqIO.SeqRecord(seq=Seq(chr_seq), id=chr_name, description=description)]
    if multiline:
        SeqIO.write(new_fasta, file_path, 'fasta')
    else:
        # output entire DNA sequence in a single line
        new_fasta = [SeqIO.FastaIO.as_fasta_2line(record) for record in new_fasta]
        with open(file_path, 'w') as handle:
            handle.writelines(new_fasta)

    return file_path


@typechecked
def save_genome_to_fasta(output_path: Path, genome: Dict[str, str], description: str = '', multiline: bool = True):
    if not output_path.exists():
        output_path.mkdir(parents=True)
    for chr_name, chr_seq in genome.items():
        save_chr_to_fasta(output_path, chr_name, chr_seq, description, multiline)


def species_reference_root(cfg: DictConfig) -> Path:
    references_dir = Path(cfg.paths.ref_dir)
    species_root_dir = references_dir / cfg.species_name
    return species_root_dir


def ref_chromosomes_path(cfg: DictConfig) -> Path:
    return species_reference_root(cfg) / 'chromosomes'


def run(cfg):
    species_path = species_reference_root(cfg)

    if species_path.exists() and not cfg.reference.overwrite:
        print(f'Reference genome for species `{cfg.species_name}` already exists. Skipping.')
        return False

    genome = get_random_genome(dict(cfg.reference.chromosomes), cfg.reference.gc_content, cfg.seed)
    save_genome_to_fasta(ref_chromosomes_path(cfg), genome, multiline=False)

    return True



@hydra.main(version_base="1.2", config_path="../../config")
def main(cfg):
    print("Running random species generator step...")
    run(cfg.reference)


if __name__ == '__main__':
    main()
