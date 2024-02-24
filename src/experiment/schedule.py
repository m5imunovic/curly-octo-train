import concurrent.futures
import logging
import subprocess
import tempfile
from pathlib import Path

from omegaconf import DictConfig
from typeguard import typechecked

import experiment.experiment_utils as eu
import reference.reference_utils as ru
from asm.assembler import run as assembly_task
from experiment.scenario_schema import Scenario, collect_all_species_defs, load_scenario
from graph.db_graph import run as graph_task
from reads.simulate_reads import run as sequencing_task
from reference.genome_generator import ensure_references_exist

logger = logging.Logger(__name__)


@typechecked
def create_sequencing_jobs(scenario: Scenario, staging_root: Path, reference_paths: dict) -> list:
    jobs = []
    job_nr = 0
    for item in scenario.items:
        # get path to the species
        # TODO: maybe batch is better name than sample
        for sample in item.samples:
            # create path to reference files and check they exist
            reference_root = reference_paths[item.species_name]
            genome = eu.get_genome_paths(ru.ref_chromosomes_path(reference_root), sample.chromosomes)
            # get the seeds and create the list of jobs
            sequencing_seeds = eu.get_sequencing_seeds(scenario.subset, sample.count, sample.init_seed)
            # chr_str = "_".join([chr.name for chr in sample.chromosomes])
            jobs.extend(
                [
                    {
                        "genome": genome,
                        "seed": seed,
                        "output_path": Path(staging_root / f"{job_nr + offset}" / "reads"),
                        # "name": f"{item.species_name.removesuffix('.yaml')}_S{seed:05d}_{chr_str}",
                    }
                    for offset, seed in enumerate(sequencing_seeds)
                ]
            )
            job_nr += len(sequencing_seeds)

    return jobs


def create_assembly_jobs(read_jobs: list) -> list:
    jobs = []
    for read_job in read_jobs:
        # find fastq files in reads
        jobs.extend(
            [
                {
                    "genome": read_job["genome"],
                    # TODO: should be dynamic fastq name based on the reads config or (better) returned from the reads job
                    "reads": [read_job["output_path"] / "sim_0001.fastq"],
                    "output_path": read_job["output_path"].parent / "assemblies",
                    "threads": 15,
                }
            ]
        )

    return jobs


@typechecked
def create_graph_jobs(assembly_jobs: list) -> list:
    jobs = []
    for assembly_job in assembly_jobs:
        # find assembly files in assemblies
        jobs.extend(
            [
                {
                    "assembly_path": assembly_job["output_path"],
                    "output_path": assembly_job["output_path"].parent / "graph",
                }
            ]
        )
    return jobs


def run_sequencing_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        cfg.reads.params.long.seed = job["seed"]
        produced_files = sequencing_task(cfg, genome=job["genome"], output_path=job["output_path"])

    return produced_files


def run_assembly_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        # Warning, this modifies the input jobs
        threads = job.pop("threads")
        cfg.asm.params.long.threads = threads
        produced_files = assembly_task(cfg, **job)
        # restore original job
        job["threads"] = threads

    return produced_files


def run_graph_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        produced_files = graph_task(cfg.graph, **job)

    return produced_files


def run_collect_job(dataset_path: Path, graph_files: dict):
    for graph_file in graph_files["artifacts"]:
        # TODO: make this more robust as we might have multiple .pt files for some reason
        if Path(graph_file).suffix == ".pt":
            idx = len(list(dataset_path.glob("*.pt")))
            src = str(graph_file)
            dst = dataset_path / f"{idx}.pt"
            logger.info(f"Copying {src} to {dst}")
            # shutil.move fails on alternative setup with SMB drive (permissions issue)
            subprocess.run(f"cp {src} {dst}", shell=True)
            break


def run_cleanup_job(step_files: dict):
    for job, file_group in step_files.items():
        for group, files in file_group.items():
            logger.info(f"Cleaning up {job} {group} files")
            for f in files:
                if f.exists():
                    logger.info(f"Removing {f}")
                    f.unlink(missing_ok=True)
                else:
                    logger.warning(f"File {f} not found")


def run(cfg: DictConfig):
    scenario = load_scenario(cfg.experiment.scenario.name)
    # get the species name and start the reference.py
    species_defs = collect_all_species_defs(scenario)
    reference_paths = ensure_references_exist(cfg, species_defs)
    if reference_paths is None:
        logger.error("Reference generation failed.")
        return

    experiment_root = eu.get_create_experiment_root(
        Path(cfg.paths.exp_dir), scenario.dataset, scenario.subset, cfg.experiment.experiment_id
    )

    dataset_path = Path(cfg.paths.datasets_dir) / scenario.dataset / scenario.subset / "raw"
    dataset_path.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as staging_dir:
        staging_dir = Path(staging_dir)
        sequencing_jobs = create_sequencing_jobs(scenario, staging_dir, reference_paths)
        assembly_jobs = create_assembly_jobs(sequencing_jobs)
        raw_dir = experiment_root.parent / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        graph_jobs = create_graph_jobs(assembly_jobs)

        def create_sample(sequencing_job, assembly_job, graph_job):
            step_files = {}
            step_files["reads"] = run_sequencing_jobs(cfg, [sequencing_job])
            step_files["assemblies"] = run_assembly_jobs(cfg, [assembly_job])
            step_files["graph"] = run_graph_jobs(cfg, [graph_job])
            return step_files

        def collect_sample(dataset_path, step_files):
            run_collect_job(dataset_path, step_files["graph"])
            run_cleanup_job(step_files)

        logger.info(f"Running {len(sequencing_jobs)} jobs in total")
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future_produced_files = {
                executor.submit(create_sample, s_job, a_job, g_job)
                for s_job, a_job, g_job in zip(sequencing_jobs, assembly_jobs, graph_jobs)
            }
            job_cnt = 0
            for future in concurrent.futures.as_completed(future_produced_files):
                step_files = future.result()
                collect_sample(dataset_path, step_files)
                logger.info(f"Finished job {job_cnt}")
                job_cnt += 1

        logger.info("Experiment finished.")
