import logging
import random
from pathlib import Path

import hydra
from Bio import SeqIO
from typeguard import typechecked

import utils.path_helpers as ph

logger = logging.getLogger(__name__)


@typechecked
def sample_reads_with_probability(reads_file: Path, sampled_file: Path, probability: float, seed: int):
    random.seed(seed)
    sampled_dir = sampled_file.parent
    if not sampled_dir.exists():
        sampled_dir.mkdir(parents=True)

    processed = 0
    selected = 0
    with open(sampled_file, "w") as out_handle:
        for record in SeqIO.parse(reads_file, "fastq"):
            processed += 1
            if processed % 10000 == 0:
                logger.info(f"Processed {processed} reads, selected {selected} samples...")
            if random.random() <= probability:
                selected += 1
                SeqIO.write(record, out_handle, "fastq")


@hydra.main(version_base=None, config_path=str(ph.get_config_root() / "reads"), config_name="sample_reads_fastq.yaml")
def main(cfg):
    sample_reads_with_probability(Path(cfg.reads_file), Path(cfg.sampled_file), cfg.probability, cfg.seed)


if __name__ == "__main__":
    main()
