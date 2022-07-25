from pathlib import Path
from typing import Dict, Optional, Sequence

import hydra
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
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


def run(cfg):
    # For now use hardcoded root path for storing genome files
    output_path = ph.get_ref_path()
    species_path = output_path / cfg.name

    if species_path.exists() and not cfg.overwrite:
        print(f'Reference genome for {cfg.name} already exists. Skipping.')
        return

    genome = get_random_genome(dict(cfg.chromosomes), cfg.gc_content, cfg.seed)
    save_genome_to_fasta(species_path / 'chromosomes', genome, multiline=False)



@hydra.main(version_base=None, config_path="../../config/reference", config_name="random_species")
def main(cfg):
    print("Running random species generator step...")
    run(cfg)


if __name__ == '__main__':
    main()
