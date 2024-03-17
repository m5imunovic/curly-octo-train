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
    probability: str | None
    chromosomes: list["Chromosome"]


@dataclass
class Chromosome:
    """Chromosome dataclass."""

    name: str


def load_scenario(scenario_name: str) -> Scenario:
    # user can supply name both with and without suffix (compatibility with yaml in hydra)
    scenario_with_suffix = f"{scenario_name.removesuffix('.json')}.json"
    scenario_path = eu.get_scenario_root() / scenario_with_suffix
    with open(scenario_path) as f:
        scenario = Scenario.from_json(f.read())

    return scenario


def collect_all_species_defs(scenario: Scenario) -> set[str]:
    references = set()
    for sample in scenario.items:
        references.add(sample.species_name)
    return references
