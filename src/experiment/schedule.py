import concurrent.futures
import fnmatch
import json
import logging
import os
import subprocess
import tempfile
import yaml
from itertools import product
from pathlib import Path

import pandas as pd
from omegaconf import DictConfig
from typeguard import typechecked

import experiment.experiment_utils as eu
import reference.reference_utils as ru
from asm.assembler import run as assembly_task
from experiment.scenario_schema import Scenario, collect_all_species_defs, load_scenario
from experiment.dataset_summary import create_entries_summary, entries_summary_to_csv
from graph.db_graph import run as graph_task
from reads.rsimulator import READ_FILE
from reads.simulate_reads import run as sequencing_task
from reference.genome_generator import ensure_references_exist

logger = logging.getLogger(__name__)


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
            reads = eu.get_read_paths(ru.ref_real_reads_path(reference_root), sample.chromosomes)
            probabilities = eu.get_sequencing_probabilities(sample.probability)
            # get the seeds and create the list of jobs
            sequencing_seeds = eu.get_sequencing_seeds(scenario.subset, sample.count, sample.init_seed)
            for probability, seed in product(probabilities, sequencing_seeds):
                jobs.append(
                    {
                        "genome": genome,
                        "seed": seed,
                        "reads": reads,
                        "output_path": Path(staging_root / f"{job_nr}" / "reads"),
                        "probability": probability,
                    }
                )
                job_nr += 1

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
                    "reads": [read_job["output_path"] / READ_FILE],
                    "output_path": read_job["output_path"].parent / "assemblies",
                }
            ]
        )

    return jobs


@typechecked
def create_graph_jobs(assembly_jobs: list, subset: str) -> list:
    jobs = []
    for assembly_job in assembly_jobs:
        # find assembly files in assemblies
        jobs.extend(
            [
                {
                    "assembly_path": assembly_job["output_path"],
                    "output_path": assembly_job["output_path"].parent / "graph",
                    "subset": subset,
                }
            ]
        )
    return jobs


def run_sequencing_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        cfg.reads.params.long.seed = job["seed"]
        produced_files = sequencing_task(cfg, **job)

    return produced_files


def run_assembly_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        produced_files = assembly_task(cfg, **job)

    return produced_files


def run_graph_jobs(cfg: DictConfig, jobs: list) -> dict:
    produced_files = {}
    for job in jobs:
        produced_files = graph_task(cfg.graph, **job)

    return produced_files


def run_collect_job(raw_path: Path, graph_files: dict) -> tuple[int, dict] | None:
    for graph_file in graph_files["artifacts"]:
        # TODO: make this more robust as we might have multiple .pt files for some reason
        if Path(graph_file).suffix == ".pt":
            idx = len(list(raw_path.glob("*.pt")))
            src = str(graph_file)
            with open(Path(graph_file).with_suffix(".json")) as fh:
                features = json.load(fh)
            dst = raw_path / f"{idx}.pt"
            logger.info(f"Copying {src} to {dst}")
            # shutil.move fails on alternative setup with SMB drive (permissions issue)
            subprocess.run(f"cp {src} {dst}", shell=True)
            return idx, features

    return None


def should_keep_file(job: str, fpath: str, keep: DictConfig | None) -> bool:
    if not keep:
        return False

    if job in keep:
        keep_patterns = keep[job]
        if any(fnmatch.fnmatch(fpath, pattern) for pattern in keep_patterns):
            return True

        if any(Path(fpath).name == pattern for pattern in keep_patterns):
            return True

    return False


def run_keep_job(dataset_root: Path, step_files: dict, idx: int, keep: DictConfig | None):
    """Runs the keep job of the files in experiment. Keeps some of the files for the evaluation and debugging purposes
    by moving them into specified directory.

    Args:
        dataset_root: root path used for saving
        step_files (dict): contains a mapping from a jobs (reads, assemblies, graph) to files.
            Files are kept as dictionaries grouped as job metadata and artifacts
        idx: index of data sample in the test dataset
        keep (DictConfig): Filtering configuration, contains artifact_path as a location go keep
        the files, either relative to dataset's train-val-test folders or absolute path location.
        and job (reads, assemblies, graph) filters as lists of patterns matchable by regex expression
    """

    for job, file_group in step_files.items():
        for group, files in file_group.items():
            if group != "artifacts":
                continue
            for f in files:
                if should_keep_file(job, str(f), keep):
                    if not f.exists():
                        logger.error(f"File {f} is missing but should be kept")
                        continue

                    relative_path = str(f).split(os.sep)[4:]  # ""/"tmp"/"tmpxyz"/ "job_id" /"relative_path..."
                    # job_id and idx might be different if we are updating existing dataset
                    output_path = (dataset_root / keep.artifacts_path / str(idx)).joinpath(*relative_path)
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    logger.info(f"Copy file {f} to {output_path}.")
                    subprocess.run(f"cp {f} {output_path}", shell=True)
                    logger.info(f"Copied file {f} to {output_path}.")


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


def save_run_metadata(cfg: DictConfig, scenario: Scenario, collected_samples: dict, dataset_root: Path) -> Path:
    metadata_path = dataset_root / "metadata"
    metadata_path.mkdir(parents=True, exist_ok=True)
    run_path = metadata_path / "run" / cfg.experiment.experiment_id
    run_path.mkdir(parents=True)
    logger.info("Saving dataset metadata")
    with open(run_path / "scenario.yaml", "w") as handle:
        yaml.safe_dump(scenario.to_dict(), handle)

    entries_summary = create_entries_summary(cfg, scenario, collected_samples)
    entries_summary_to_csv(entries_summary, run_path)
    return run_path


def merge_run_metadata(subset_path: Path, run_path: Path):
    new_summary = pd.read_csv(run_path / "summary.csv", index_col=False)
    subset_summary_path = subset_path / "summary.csv"
    if subset_summary_path.exists():
        subset_summary = pd.read_csv(subset_summary_path, index_col=False)
        latest_revision = subset_summary["Revision"].astype("int").max()
        new_summary["Revision"] = latest_revision + 1
        new_summary = subset_summary.append(new_summary, ignore_index=True)
    else:
        new_summary["Revision"] = 1

    new_summary.to_csv(subset_summary_path, index=False, float_format="%.3f")


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

    dataset_root = Path(cfg.paths.datasets_dir) / scenario.dataset
    dataset_root.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as staging_dir:
        raw_path = dataset_root / scenario.subset / "raw"
        raw_path.mkdir(parents=True, exist_ok=True)
        staging_dir = Path(staging_dir)
        sequencing_jobs = create_sequencing_jobs(scenario, staging_dir, reference_paths)
        assembly_jobs = create_assembly_jobs(sequencing_jobs)
        raw_dir = experiment_root.parent / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        graph_jobs = create_graph_jobs(assembly_jobs, scenario.subset)

        def create_sample(sequencing_job, assembly_job, graph_job):
            step_files = {}
            step_files["reads"] = run_sequencing_jobs(cfg, [sequencing_job])
            step_files["assemblies"] = run_assembly_jobs(cfg, [assembly_job])
            step_files["graph"] = run_graph_jobs(cfg, [graph_job])
            return step_files

        def collect_sample(raw_path, step_files, cfg):
            idx, features = run_collect_job(raw_path, step_files["graph"])
            if scenario.subset == "test":
                run_keep_job(raw_path.parent.parent, step_files, idx, cfg.experiment.keep)
            run_cleanup_job(step_files)
            return idx, features

        logger.info(f"Running {len(sequencing_jobs)} jobs in total")
        collected_samples = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=cfg.experiment.max_workers) as executor:
            future_produced_files = {
                executor.submit(create_sample, s_job, a_job, g_job)
                for s_job, a_job, g_job in zip(sequencing_jobs, assembly_jobs, graph_jobs)
            }
            job_cnt = 0
            for future in concurrent.futures.as_completed(future_produced_files):
                step_files = future.result()
                idx, features = collect_sample(raw_path, step_files, cfg)
                collected_samples[idx] = features
                logger.info(f"Finished job {job_cnt}")
                job_cnt += 1

    run_path = save_run_metadata(cfg, scenario, collected_samples, dataset_root)
    merge_run_metadata(dataset_root / scenario.subset, run_path)
    experiment_root.rmdir()
    logger.info("Experiment finished.")
