import logging
import random
import subprocess
from pathlib import Path

import hydra
from Bio import SeqIO
from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from reads.reads_simulator import simulator_factory
from utils.io_utils import compose_cmd_params

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


def run(cfg: DictConfig):
    assert Path(cfg.reads.params.long.genome).exists(), f"Genome file {cfg.reads.params.long.genome} does not exist!"
    assert str(cfg.sample.seed).isdigit(), f"Seed {cfg.sample.seed} is not a number!"
    # TODO: assert sample.origin exists and add printout
    assert Path(cfg.sample.reads_file).exists(), f"Reads file {cfg.sample.reads_file} does not exist!"

    sample_reads_with_probability(
        Path(cfg.sample.reads_file), Path(cfg.sample.sampled_file), cfg.sample.probability, int(cfg.sample.seed)
    )
    # TODO: refactor PbSim3 to be able to create profile as config option
    pbsim3 = simulator_factory("pbsim3", cfg.reads)
    option_params = compose_cmd_params(cfg.reads.params)
    cmds = []
    cmds.append(f"{pbsim3.simulator_exec} {option_params}")
    cmds.append("rm sd_0001.maf")
    cmds.append("rm sd_0001.fastq")
    cmds.append("rm sd_0001.ref")
    # TODO: remove the intermediate fastq file
    for cmd in cmds:
        logger.info(f"RUN::profile:: {cmd}")
        subprocess.run(cmd, shell=True, cwd=cfg.reads.profile.path)
    # TODO: add wandb upload


@hydra.main(version_base=None, config_path=str(ph.get_config_root() / "reads"), config_name="sample_reads_fastq.yaml")
def main(cfg):
    run(cfg)


if __name__ == "__main__":
    main()
