from typing import Sequence

import numpy as np
from typeguard import typechecked


@typechecked
def get_random_chromosome(
    sequence_length: int, rng: np.random.Generator, base_p: Sequence[float] | None = None
) -> str:
    bases = rng.choice(["A", "C", "G", "T"], size=sequence_length, p=base_p)
    return "".join(bases)


@typechecked
def get_random_reference(chromosomes: dict[str, int], seed: int, gc_content: float | None = None) -> dict[str, str]:
    rng = np.random.default_rng(12345)
    np.random.seed(seed)
    base_p = [(1 - gc_content) / 2, gc_content / 2, gc_content / 2, (1 - gc_content) / 2] if gc_content else None

    genome = {}
    for chr_name, chr_length in chromosomes.items():
        chr_seq = get_random_chromosome(chr_length, rng, base_p)
        genome[chr_name] = chr_seq

    return genome
