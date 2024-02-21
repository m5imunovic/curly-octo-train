from omegaconf import DictConfig
from typeguard import typechecked

from reads.reads_simulator import simulator_factory
from utils.io_utils import get_job_outputs


@typechecked
def run(cfg: DictConfig, **kwargs) -> dict:
    # species_filter contains species name (unique path identifier in the reference directory)
    # chromosomes from which to sequence, how many reads to simulate and range of nucleotides

    exec_args = {
        # Top level output path
        "output_path": None,
        # Path to the reference genome directory (can contain one or multiple files)
        "genome": None,
    }
    exec_args.update(kwargs)
    simulator = simulator_factory(simulator=cfg.reads.name, cfg=cfg.reads)

    simulator(**exec_args)

    return get_job_outputs(exec_args)
