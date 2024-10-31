import numpy as np


BASES = ["A", "C", "T", "G"]


def random_dna_sequence(length):
    return "".join(np.random.choice(BASES, size=length).tolist())


def insert_element(sequence, element, frequency):
    """Insert the element into the sequence at random positions at a given frequency."""
    num_insertions = int(len(sequence) * frequency)

    positions = list(sorted(np.random.choice(len(sequence), num_insertions, replace=False)))
    positions.append(len(sequence))

    parts = [sequence[: positions[0]]]
    for pos_start, pos_end in zip(positions[:-1], positions[1:]):
        parts.append(element)
        parts.append(sequence[pos_start:pos_end])

    return "".join(parts)


def generate_repeat_element(length):
    """
    Function to generate a random repeat element (e.g. LINEs, SINEs)
    - LINEs are anywhere between 0.5 and 8 kB long but more typically about 1kB,
    and there are around 100 of those that are around 6 kB in human chromosome.
    - SINEs elements are 100B to 300B long.
    - Segmental Duplications are between 1 and 400 kB length
    """
    return random_dna_sequence(length)


def generate_telomeric_sequence(repeat_unit="TTAGGG", repeat_count=500):
    """Simulate telomeric sequence."""
    return repeat_unit * repeat_count


def generate_centromeric_repeat(length=50000):
    """Simulate centromeric repeat (e.g., alpha-satellite)"""
    return "ACGT" * (length // 4)  # Simulate alpha-satellite repeat


def generate_microsatellite(repeat_unit="CAT", repeat_length=20):
    """Simulate microsatellites (short tandem repeats).

    These are 1 to 6 bases long and repeat typically 5 to 50 times
    """
    return repeat_unit * (repeat_length // len(repeat_unit))


def generate_repeat_sequence_complex(chromosome):
    """Generate sequence with repeats."""

    print("Adding centromere")
    centromeric_seq = generate_centromeric_repeat()
    middle = len(chromosome) // 2
    chromosome = chromosome[:middle] + centromeric_seq + chromosome[middle:]

    print("Adding LINEs")
    line_seq = generate_repeat_element(3000)
    chromosome = insert_element(chromosome, line_seq, frequency=0.001)

    print("Adding SINEs")
    sine_seq = generate_repeat_element(200)
    chromosome = insert_element(chromosome, sine_seq, frequency=0.002)

    print("Generating microsatelites")
    microsatelite_seq = generate_microsatellite()
    chromosome = insert_element(chromosome, microsatelite_seq, frequency=0.005)

    print("Generating segmental duplication")
    segmental_duplication = generate_repeat_element(3000)
    chromosome = insert_element(chromosome, segmental_duplication, frequency=0.0001)

    print("Generating telomeres")
    telometric_seq = generate_telomeric_sequence()
    chromosome = telometric_seq + chromosome + telometric_seq

    print(f"Total length of a chromosome is {len(chromosome)}")

    return chromosome
