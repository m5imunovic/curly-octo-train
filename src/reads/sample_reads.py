import logging
import random
from datetime import datetime
from pathlib import Path

from Bio import SeqIO
from typeguard import typechecked

from utils.io_utils import ensure_parent_dir_exists

logger = logging.getLogger(__name__)


@typechecked
def sample_reads_with_probability(reads_file: Path, sampled_file: Path, probability: float, seed: int) -> int:
    """Sample reads from a FASTQ file based on a given probability. The purpose is to reduce input reads file before
    using it to create a simulation profile.

    Args:
        reads_file (Path): Path to the input FASTQ file.
        sampled_file (Path): Path to the output file where the sampled reads will be written.
        probability (float): Probability of selecting each read.
        seed (int): Seed value for the random number generator.

    Returns:
        int: Number of selected reads.
    """
    ensure_parent_dir_exists(sampled_file)

    start = datetime.now()
    with open(sampled_file, "a") as out_handle:
        random.seed(seed)
        processed, selected = 0, 0
        for record in SeqIO.parse(reads_file, "fastq"):
            processed += 1
            if processed % 10000 == 0:
                logger.info(f"Processed {processed} reads, selected {selected} samples...")
            if random.random() <= probability:
                selected += 1
                SeqIO.write(record, out_handle, "fastq")

        logger.info(f"Finished sampling {processed} reads from {reads_file} after {datetime.now() - start}!")
        logger.info(f"Selected {selected} samples")

    return selected
