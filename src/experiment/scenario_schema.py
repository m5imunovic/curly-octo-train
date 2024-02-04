"""
Generated with dataclass-wizard -> cat scenario.json | wiz gs - scenario_schema
"""
from dataclasses import dataclass

from dataclass_wizard import JSONWizard

import experiment.experiment_utils as eu


@dataclass
class Scenario(JSONWizard):
    """Data dataclass."""

    dataset: str
    subset: str
    items: list["Item"]


@dataclass
class Item:
    """Item dataclass."""

    species_name: str
    samples: list["Sample"]


@dataclass
class Sample:
    """Sample dataclass."""

    init_seed: int
    count: int
    chromosomes: list["Chromosome"]


@dataclass
class Chromosome:
    """Chromosome dataclass."""

    name: str


def load_scenario(scenario_name: str) -> Scenario:
    scenario_path = eu.get_scenario_root() / f"{scenario_name}.json"
    with open(scenario_path) as f:
        scenario = Scenario.from_json(f.read())

    return scenario


def collect_all_species_defs(scenario: Scenario) -> set[str]:
    references = set()
    for sample in scenario.items:
        references.add(sample.species_name)
    return references
