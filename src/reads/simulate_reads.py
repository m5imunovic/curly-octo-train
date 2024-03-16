from pathlib import Path
from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from reads.rsimulator import RSimulator
from reads.reads_simulator import PbSim3
from reads.sampler import ReadSampler
from utils.io_utils import get_job_outputs


@typechecked
def simulator_factory(simulator: str, cfg: DictConfig) -> RSimulator:
    vendor_dir: Path = ph.get_vendor_path() if "paths" in cfg else ph.get_vendor_path()
    if simulator == "pbsim3":
        return PbSim3(cfg=cfg, vendor_dir=vendor_dir)
    if simulator == "sampler":
        return ReadSampler(cfg, vendor_dir=None)
    raise ValueError(f"Unknown simulator name {simulator}")


@typechecked
def run(cfg: DictConfig, **kwargs) -> dict:

    exec_args = {
        # Top level output path
        "output_path": None,
        # Path to the reference genome directory (can contain one or multiple files)
        "genome": None,
        "probability": 1.0,
        "reads": None,
    }
    exec_args.update(kwargs)

    simulator = simulator_factory(simulator=cfg.reads.name, cfg=cfg.reads)
    simulator(**exec_args)

    return get_job_outputs(exec_args)
