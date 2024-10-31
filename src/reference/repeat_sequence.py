from pathlib import Path

import numpy as np
from Bio import SeqIO

from reference.repeat_sequence_complex import generate_repeat_sequence_complex


def select_random_sequence(fasta_file, seq_length=1_000_000):
    """Function to select random sequence of ~1Mb from the reference FASTA file."""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    genome_length = len(record.seq)
    start = np.random.randint(0, genome_length - seq_length)
    end = start + seq_length

    return record.seq[start:end]


def add_repeats(sequence, chunk_low=5000, chunk_hi=6000, range_low=1000, range_hi=3000, rep_low=1, rep_hi=5):
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
        rep_cnt = np.random.randint(rep_low, rep_hi)
        for rep in range(rep_cnt):
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


def get_repeat_sequence(cfg, output_path):
    fasta_file = Path(cfg.chr_src_path)
    assert fasta_file.exists(), f"{fasta_file} does not exists"
    seed = cfg.chr_src_seed
    np.random.seed(seed)
    output_path = Path(output_path)
    chromosome_path = output_path / "chromosomes"
    chromosome_path.mkdir(parents=True, exist_ok=True)

    mode = cfg.mode
    sequence_len = cfg.sequence_len
    mr = cfg.mutation_rates
    mutation_rates = np.arange(mr.start, mr.stop, mr.step)
    p_choices = cfg.p_choices
    # Step 1: Select a random sequence of ~1Mb
    for idx, mutation_rate in enumerate(mutation_rates):
        print(f"Generating {idx+1} of {len(mutation_rates)}")
        sequence = select_random_sequence(fasta_file, sequence_len)

        # Step 2: Add repeats to the sequence
        if mode == "simple":
            modified_sequence = add_repeats(sequence)
        elif mode == "complex":
            modified_sequence = generate_repeat_sequence_complex(str(sequence))

        # Step 3: Generate two haplotypes with random mutations
        print("Generating haplotypes...")
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
