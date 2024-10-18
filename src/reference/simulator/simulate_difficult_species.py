import json
from pathlib import Path

import numpy as np
from Bio import SeqIO


def select_random_sequence(fasta_file, seq_length=1_000_000):
    """Function to select random sequence of ~1Mb from the reference FASTA file."""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    genome_length = len(record.seq)
    start = np.random.randint(0, genome_length - seq_length)
    end = start + seq_length

    return record.seq[start:end]


def add_repeats(sequence, chunk_low=5000, chunk_hi=6000, range_low=1000, range_hi=3000):
    """Introduce repeats of 1-3Kb every 5-6Kb."""
    new_sequence = []
    start = 0
    end = len(sequence)
    while start < end:
        # Handle last chunk
        if start + chunk_hi >= end:
            new_sequence.append(str(sequence[start:]))
            break

        # Extract the next chunk (5-6Kb)
        chunk_size = np.random.randint(chunk_low, chunk_hi)
        chunk = sequence[start : start + chunk_size]
        new_sequence.append(str(chunk))

        # Add a random repeat of the last 1-3Kb from the current chunk
        repeat_length = np.random.randint(range_low, range_hi)
        repeat = chunk[-repeat_length:]
        new_sequence.append(str(repeat))

        start += chunk_size

    return "".join(new_sequence)


def introduce_mutations(sequence, mutation_rate, p_choices):
    """Introduce random mutations with a given rate (0.1-1%)"""
    sequence = list(sequence)
    num_mutations = int(len(sequence) * mutation_rate)

    # delete, insert, mutate
    choices = ("D", "I", "M")
    mutations = {}
    for _ in range(num_mutations):
        # Choose a random position
        pos = np.random.randint(0, len(sequence) - 1)
        choice = np.random.choice(choices, p=p_choices)
        mutations[pos] = choice

    mutations = dict(sorted(mutations.items()))
    cigar = "".join([f"{pos}{mut}" for pos, mut in mutations.items()])
    mutated = []
    for i in range(len(sequence)):
        if i in mutations:
            mutation = mutations[i]
            if mutation == "D":
                continue
            # either mutate to other base or insert any
            bases = ["A", "T", "C", "G"]
            if mutation == "M":
                bases.remove(sequence[pos])  # Remove the original base
            new_base = np.random.choice(bases)
            mutated.append(new_base)
        else:
            # keep original
            mutated.append(sequence[i])

    return "".join(mutated), cigar


def generate_haplotypes(sequence, mutation_rate, p_choices):
    "Generate two haplotypes with mutations"
    haplotype1 = introduce_mutations(sequence, mutation_rate, p_choices)
    haplotype2 = introduce_mutations(sequence, mutation_rate, p_choices)
    return haplotype1, haplotype2


def run(cfg):
    fasta_file = Path(cfg.chromosome_path)
    assert fasta_file.exists(), f"{fasta_file} does not exists"
    seed = cfg.seed
    np.random.seed(seed)
    output_path = Path(cfg.output_path)
    chromosome_path = output_path / "chromosomes"
    chromosome_path.mkdir(parents=True, exist_ok=True)

    with open(output_path / "species_info.json", "w") as f:
        json.dump(cfg, f, indent=4)

    sequence_len = cfg.sequence_len
    mutation_rates = cfg.mutation_rates
    p_choices = cfg.p_choices
    # Step 1: Select a random sequence of ~1Mb
    for idx, mutation_rate in enumerate(mutation_rates):
        sequence = select_random_sequence(fasta_file, sequence_len)

        # Step 2: Add repeats to the sequence
        modified_sequence = add_repeats(sequence)

        # Step 3: Generate two haplotypes with random mutations
        haplotype1, haplotype2 = generate_haplotypes(modified_sequence, mutation_rate, p_choices)
        name_mat = f"chr{idx+1}_MATERNAL"
        with open((chromosome_path / name_mat).with_suffix(".fasta"), "w") as f:
            f.write(f">{name_mat}\n")
            f.write(haplotype1[0])
            f.write("\n")
        name_pat = f"chr{idx+1}_PATERNAL"
        with open((chromosome_path / name_pat).with_suffix(".fasta"), "w") as f:
            f.write(f">{name_pat}\n")
            f.write(haplotype2[0])
            f.write("\n")


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self


if __name__ == "__main__":
    cfg = {
        "chromosome_path": "/data/references/homo_sapiens/chm13_v2/chromosomes/chr10.fasta",
        "mutation_rates": np.arange(0.001, 0.0105, 0.001).tolist(),
        "seed": 400,
        "sequence_len": 1000000,
        "p_choices": (0.5, 0.25, 0.25),
        "output_path": "/data/references/repeat_species/v1/",
    }
    run(AttrDict(cfg))
