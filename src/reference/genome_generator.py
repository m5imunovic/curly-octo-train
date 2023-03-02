import json
import random
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Sequence

import hydra
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from omegaconf import DictConfig
from typeguard import typechecked

rng = np.random.default_rng(12345)


@typechecked
def get_random_chromosome(
    sequence_length: int, base_p: Optional[Sequence[float]] = None, seed: Optional[int] = None
) -> str:
    if seed:
        np.random.seed(seed)
    bases = rng.choice(["A", "C", "G", "T"], size=sequence_length, p=base_p)
    return "".join(bases)


@typechecked
def get_random_genome(
    chromosomes: Dict[str, int], gc_content: Optional[float] = None, seed: Optional[int] = None
) -> Dict[str, str]:

    base_p = [(1 - gc_content) / 2, gc_content / 2, gc_content / 2, (1 - gc_content) / 2] if gc_content else None

    genome = {}
    seed_chr_offset = 0
    if seed is None:
        seed = random.randint(0, 1000000)
        print(f"Using randomly generated {seed=}.")

    for chr_name, chr_length in chromosomes.items():
        chr_seq = get_random_chromosome(chr_length, base_p, seed + seed_chr_offset)

        genome[chr_name] = chr_seq
        seed_chr_offset += 1

    return genome


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
def save_genome_to_fasta(output_path: Path, genome: Dict[str, str], description: str = "", multiline: bool = True):
    if not output_path.exists():
        output_path.mkdir(parents=True)
    for chr_name, chr_seq in genome.items():
        save_chr_to_fasta(output_path, chr_name, chr_seq, description, multiline)


def species_reference_root(cfg: DictConfig) -> Path:
    references_dir = Path(cfg.paths.ref_dir)
    species_name = str(cfg.species_name)
    return references_dir / species_name


def ref_chromosomes_path(cfg: DictConfig) -> Path:
    return species_reference_root(cfg) / "chromosomes"


def run(cfg):
    species_path = species_reference_root(cfg)

    if species_path.exists():
        print(f"Reference genome for species `{cfg.species_name}` already exists at location: \n{species_path}")
        print("Skipping...")
        return False

    try:
        print(cfg)
        species_path.mkdir(parents=True)

        species_info = {
            "species_name": cfg.species_name,
            "chromosomes": dict(cfg.reference.chromosomes),
            "gc_content": cfg.reference.gc_content,
            "seed": cfg.seed,
            "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        with open(species_path / "species_info.json", "w") as handle:
            json.dump(species_info, handle, indent=4)

            genome = get_random_genome(dict(cfg.reference.chromosomes), cfg.reference.gc_content, cfg.seed)
            save_genome_to_fasta(ref_chromosomes_path(cfg), genome, multiline=False)

    except Exception as e:
        shutil.rmtree(species_path)
        print(f"Error while generating reference genome for species `{cfg.species_name}`: {e}")
        return False

    return True


@hydra.main(version_base=None, config_path="../../config")
def main(cfg):
    # TODO: use logger instead
    print("Starting random species generator step...")
    run(cfg.reference)
    print("Finished species generator step.")


if __name__ == "__main__":
    main()
