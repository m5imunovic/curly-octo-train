from itertools import product
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from omegaconf import DictConfig

import experiment.experiment_utils as eu
from experiment.scenario_schema import Scenario, SCHEMA


@dataclass
class DatasetInfoEntry:
    Subset: str
    Id: int
    Species: str
    Chromosomes: list[str]
    Nodes: int
    Edges: int
    Seed: int
    Probability: float
    ProfileId: str | None
    Depth: str | None
    RunId: str
    Schema: str


def create_entries_summary(cfg: DictConfig, scenario: Scenario, collected_samples: list) -> list:
    # first constants
    subset = scenario.subset
    read_sampling = "simulate" if cfg.reads.name == "pbsim3" else "sample"
    profile_id = cfg.reads.params.long["sample-profile-id"] if read_sampling == "simulate" else ""
    depth = cfg.reads.params.long["depth"] if read_sampling == "simulate" else ""
    run_id = cfg.experiment.experiment_id
    schema = SCHEMA

    summary = []
    id = 0
    for item in scenario.items:
        species = item.species_name
        for sample in item.samples:
            chromosomes = [chr.name for chr in sample.chromosomes]
            probabilities = eu.get_sequencing_probabilities(sample.probability)
            sequencing_seeds = eu.get_sequencing_seeds(scenario.subset, sample.count, sample.init_seed)
            for probability, seed in product(probabilities, sequencing_seeds):
                if id < len(collected_samples):
                    entry = collected_samples[id]
                    id += 1
                    if entry is None:
                        continue
                    sample_id, features = entry
                    num_nodes, num_edges = features["num_nodes"], features["num_edges"]
                    entry = DatasetInfoEntry(
                        subset,
                        sample_id,
                        species,
                        chromosomes,
                        num_nodes,
                        num_edges,
                        seed,
                        probability,
                        profile_id,
                        depth,
                        run_id,
                        schema,
                    )
                    summary.append(entry)

    return summary


def entries_summary_to_csv(entries: list, output_path: Path, filename: str = "summary.csv"):
    if len(entries) == 0:
        return
    data = [asdict(entry) for entry in entries]
    df = pd.DataFrame(data)
    df.to_csv(output_path / filename, index=False, float_format="%.3f")
