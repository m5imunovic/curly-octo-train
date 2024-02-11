import logging
import shutil
import tempfile
from pathlib import Path

from omegaconf import DictConfig, OmegaConf
from typeguard import typechecked

import experiment.experiment_utils as eu
import reference.reference_utils as ru
import utils.path_helpers as ph
from asm.assembler import run as assembly_task
from experiment.scenario_schema import Scenario, collect_all_species_defs, load_scenario
from graph.db_graph import run as graph_task
from reads.simulate_reads import run as sequencing_task
from reference.genome_generator import run as reference_task

logger = logging.Logger(__name__)


def ensure_references_exist(cfg: DictConfig, scenario: Scenario) -> dict | None:
    """Ensures that all species references exist for the given scenario.

    Args:
        scenario (Scenario): The scenario for which to ensure species references exist.

    Returns:
        bool: True if all species references exist, False otherwise.
    """
    species_defs = collect_all_species_defs(scenario)
    config_root = ph.get_config_root()
    try:
        eu.ensure_species_def_exists(config_root, species_defs)
    except AssertionError as e:
        logger.error(f"Species definition files error:\n {e}")
        return False

    # TODO: this is a candidate for parallelization
    species_paths = {}
    species_config_root = config_root / "reference" / "species"
    for species_def in species_defs:
        cfg.reference.species = OmegaConf.load(species_config_root / f"{species_def}")
        reference_path = reference_task(cfg)
        if not reference_path:
            return None
        species_paths[species_def] = reference_path

    return species_paths


@typechecked
def create_sequencing_jobs(scenario: Scenario, sequencing_root: Path, reference_paths: dict) -> list:
    jobs = []
    for item in scenario.items:
        # get path to the species
        # TODO: maybe batch is better name than sample
        for sample in item.samples:
            # create path to reference files and check they exist
            reference_root = reference_paths[item.species_name]
            genome = eu.get_genome_paths(ru.ref_chromosomes_path(reference_root), sample.chromosomes)
            # get the seeds and create the list of jobs
            sequencing_seeds = eu.get_sequencing_seeds(scenario.subset, sample.count, sample.init_seed)
            chr_str = "_".join([chr.name for chr in sample.chromosomes])
            jobs.extend(
                [
                    {
                        "genome": genome,
                        "seed": seed,
                        "output_path": sequencing_root
                        / f"{item.species_name.removesuffix('.yaml')}_S{seed:05d}_{chr_str}",
                    }
                    for seed in sequencing_seeds
                ]
            )

    return jobs


def create_assembly_jobs(read_jobs: list, assembly_root: Path) -> list:
    jobs = []
    for read_job in read_jobs:
        # find fastq files in reads
        jobs.extend(
            [
                {
                    "genome": read_job["genome"],
                    # TODO: should be dynamic fastq name based on the reads config or (better) returned from the reads job
                    "reads": [read_job["output_path"] / "sim_0001.fastq"],
                    "output_path": assembly_root / f"{read_job['output_path'].name}",
                    "threads": 15,
                }
            ]
        )

    return jobs


@typechecked
def create_graph_jobs(assembly_jobs: list, graph_root: Path, raw_dir: Path) -> list:
    jobs = []
    # TODO: add exception clause here
    file_stems = [int(f.stem) for f in raw_dir.glob("*.pt")]
    idx = max(file_stems) + 1 if file_stems else 0
    for assembly_job in assembly_jobs:
        # find assembly files in assemblies
        jobs.extend(
            [
                {
                    "assembly_path": assembly_job["output_path"],
                    "output_path": graph_root / f"{assembly_job['output_path'].name}",
                    "idx": idx,
                }
            ]
        )
        idx += 1
    return jobs


def run_sequencing_jobs(cfg: DictConfig, jobs: list) -> list:
    for job in jobs:
        cfg.reads.params.long.seed = job["seed"]
        sequencing_task(cfg, genome=job["genome"], output_path=job["output_path"])

    return []


def run_assembly_jobs(cfg: DictConfig, jobs: list):
    for job in jobs:
        # Warning, this modifies the input jobs
        threads = job.pop("threads")
        cfg.asm.params.long.threads = threads
        assembly_task(cfg, **job)
        # restore original job
        job["threads"] = threads

    return []


def run_graph_jobs(cfg: DictConfig, jobs: list) -> list:
    for job in jobs:
        graph_task(cfg.graph, **job)

    return []


def run_cleanup_job(read_job: dict, assembly_job: dict):
    shutil.rmtree(read_job["output_path"])
    shutil.rmtree(assembly_job["output_path"])


def run(cfg: DictConfig):
    scenario = load_scenario(cfg.experiment.scenario.name)
    # get the species name and start the reference.py
    reference_paths = ensure_references_exist(cfg, scenario)
    if reference_paths is None:
        logger.error("Reference generation failed.")
        return

    experiment_root = eu.get_create_experiment_root(
        Path(cfg.paths.exp_dir), scenario.dataset, scenario.subset, cfg.experiment.experiment_id
    )

    tempfile.TemporaryDirectory()

    sequencing_jobs = create_sequencing_jobs(scenario, experiment_root / "sequencing", reference_paths)
    assembly_jobs = create_assembly_jobs(sequencing_jobs, experiment_root / "assembly")
    raw_dir = experiment_root.parent / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    graph_jobs = create_graph_jobs(assembly_jobs, experiment_root / "graph", raw_dir)

    for i in range(len(sequencing_jobs)):
        logger.info(f"Running job {i+1} of {len(sequencing_jobs)}")
        run_sequencing_jobs(cfg, [sequencing_jobs[i]])
        run_assembly_jobs(cfg, [assembly_jobs[i]])
        run_graph_jobs(cfg, [graph_jobs[i]])
        run_cleanup_job(sequencing_jobs[i], assembly_jobs[i])
