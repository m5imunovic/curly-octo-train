from dataclasses import dataclass
from pathlib import Path

from omegaconf import DictConfig, OmegaConf

from reads.reads_simulator import simulator_factory


@dataclass
class SampleProfile:
    genome: list[Path]
    seed: int


def run(cfg: DictConfig, **kwargs):
    # species_filter contains species name (unique path identifier in the reference directory)
    # chromosomes from which to sequence, how many reads to simulate and range of nucleotides

    exec_args = {
        # Top level output path
        "simulated_reads_path": None,
        # Path to the reference genome directory (can contain one or multiple files)
        "genome": None,
    }
    exec_args.update(kwargs)
    simulator = simulator_factory(simulator=cfg.reads.name, cfg=cfg.reads)

    simulator.pre_simulation_step(**exec_args)
    simulator.run(**exec_args)
    simulator.post_simulation_step(**exec_args)

    metadata_path = exec_args["simulated_reads_path"] / "metadata"
    metadata_path.mkdir(parents=True, exist_ok=True)

    with open(metadata_path / "reads.yaml", "w") as f:
        OmegaConf.save(config=cfg, f=f, resolve=True)
