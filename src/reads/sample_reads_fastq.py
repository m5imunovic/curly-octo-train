import logging
import subprocess
from pathlib import Path

from omegaconf import DictConfig

from reads.sample_reads import sample_reads_with_probability
from reads.simulate_reads import simulator_factory
from reads.sequencing_utils import upload_pbsim3_profile_to_wandb
from reference.genome_generator import ensure_references_exist
from utils.io_utils import compose_cmd_params

logger = logging.getLogger(__name__)


def check_sample_reads_config(cfg: DictConfig):
    """Check the configuration for sample reads.

    Args:
        cfg (DictConfig): The configuration dictionary.

    Raises:
        AssertionError: If the origin is not defined in the config.
        AssertionError: If the seed is not a number.
        AssertionError: If the reads file from which to which we want to subsample does not exist.
    """
    assert cfg.origin is not None, "Origin is not defined in the config!"
    assert str(cfg.seed).isdigit(), f"Seed {cfg.seed} is not a number!"
    assert Path(cfg.reads_file).exists(), f"Reads file {cfg.reads_file} does not exist!"


def check_prerequisites(cfg: DictConfig):
    """Check the prerequisites for creating a profile. Do the checking early to avoid running the whole pipeline.
    Ensure that the reference genome exists used for simulation exists.

    Args:
        cfg (DictConfig): The configuration object.

    Raises:
        AssertionError: If sample.reads contains config errors
        AssertionError: If sample.genome file does not exist
        AssertionError: If the genome file specified in the configuration does not exist.
    """
    check_sample_reads_config(cfg.sample.reads)
    ensure_references_exist(cfg, {cfg.sample.genome.species_name})
    assert Path(cfg.reads.params.long.genome).exists(), f"Genome file {cfg.reads.params.long.genome} does not exist!"


def prepare_sample_reads(cfg: DictConfig):
    """Small wrapper to reduce cfg... typing.

    Args:
        cfg (DictConfig): Only cfg.sample.reads is used.
    """
    sample_reads_with_probability(Path(cfg.reads_file), Path(cfg.sampled_file), cfg.probability, int(cfg.seed))


def pbsim3_create_profile(cfg: DictConfig):
    # TODO: refactor PbSim3 to be able to create profile as config option
    pbsim3 = simulator_factory("pbsim3", cfg.reads)
    option_params = compose_cmd_params(cfg.reads.params)
    cmds = []
    cmds.append(f"{pbsim3.simulator_exec} {option_params}")
    cmds.append("rm sd_0001.maf")
    cmds.append("rm sd_0001.fastq")
    cmds.append("rm sd_0001.ref")
    cmds.append(f"rm {cfg.sample.reads.sampled_file}")
    for cmd in cmds:
        logger.info(f"RUN::profile:: {cmd}")
        subprocess.run(cmd, shell=True, cwd=cfg.reads.profile.path)


def run(cfg: DictConfig):
    check_prerequisites(cfg)
    prepare_sample_reads(cfg.sample.reads)
    pbsim3_create_profile(cfg)
    upload_pbsim3_profile_to_wandb(cfg)
