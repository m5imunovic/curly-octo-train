from dataclasses import dataclass
from pathlib import Path

import yaml
from omegaconf import DictConfig, OmegaConf

from reads.reads_simulator import simulator_factory


@dataclass
class SampleProfile:
    genome: list[Path]
    seed: int


def run(cfg: DictConfig, **kwargs):
    # species_filter is used to control which species are simulated
    # species_filter contains species name (unique path identifier in the reference directory)
    # chromosomes from which to sequence, how many reads to simulate and range of nucleotides

    output_path = Path(cfg.paths.reads_dir) / cfg.date_mm_dd / f"S{cfg.seed}"
    exec_args = {
        # Top level output path
        "simulated_species_path": output_path,
        # Path to the reference genome directory (can contain one or multiple fasta files)
        "genome": Path(cfg.paths.ref_dir) / cfg.species_name["name"],
    }
    exec_args.update(kwargs)
    simulator = simulator_factory(simulator=cfg.reads.name, cfg=cfg.reads)

    simulator.pre_simulation_step(**exec_args)
    simulator.run(**exec_args)
    simulator.post_simulation_step(**exec_args)

    metadata_path = output_path.parent / "metadata"
    metadata_path.mkdir(parents=True, exist_ok=True)

    with open(metadata_path / f"S{cfg.seed}.yaml", "w") as f:
        OmegaConf.resolve(cfg)
        cfg_cont = OmegaConf.to_container(cfg)
        cfg_cont.pop("paths", None)
        yaml.dump(cfg_cont, f)
