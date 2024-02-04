from pathlib import Path

import utils.path_helpers as ph


def get_scenario_root() -> Path:
    return ph.get_config_root() / "experiment" / "scenarios"


def get_create_experiment_root(exp_dir: Path, dataset: str, subset: str, exp_id: str) -> Path:
    """Get or create the root directory for an experiment. Dataset and subset info comes from experiment scenario
    definition. The experiment ID comes from the app config. It is created from date and time of the experiment start,
    in YYYYMMDD_HHMM format so typically it is unique.

    Args:
        exp_dir (Path): The root directory for all the experiments.
        dataset (str): The dataset name.
        subset (str): The subset name.
        exp_id (str): The experiment ID.

    Returns:
        Path: The path to this experiment root directory.
    """
    experiment_root = exp_dir / dataset / subset / exp_id
    experiment_root.mkdir(parents=True, exist_ok=True)
    return experiment_root


def ensure_species_def_exists(config_root: Path, species_defs: list[str]):
    """Ensures that the species definition files exist in the specified directory. It is always assumed that the
    species definition files are located in the `config_root/reference/species` directory.

    Args:
        config_root (Path): The root directory of the configuration files.
        species_defs (list[str]): A list of species definition file names with file extension.

    Raises:
        AssertionError: If any of the species definition files do not exist or does not have a .yaml extension.
    """
    for species_def in species_defs:
        assert species_def.endswith(".yaml"), f"Species definition file must have a .yaml extension: {species_def}"
        species_def_path = config_root / "reference" / "species" / species_def
        assert species_def_path.exists(), f"Species definition does not exist at {species_def_path}"


def get_genome_paths(ref_chromosomes_path: Path, chromosomes: list) -> list[Path]:
    chr_paths = []
    for chr in chromosomes:
        chr_path = ref_chromosomes_path / f"{chr.name}.fasta"
        if not chr_path.exists():
            raise FileNotFoundError(f"Chromosome file does not exist at {chr_path}")
        chr_paths.append(chr_path)

    return chr_paths


def get_sequencing_seeds(subset: str, count: int, init_seed: int) -> list[int]:
    """Get a list of seeds for the given subset and count. The seeds are generated by adding the count to the initial
    seed and offset. Offset is dependent on the subset name.

    Args:
        subset (str): The subset name.
        count (int): The number of seeds to generate.
        init_seed (int): The initial seed.

    Returns:
        list[int]: A list of seeds.
    Raises:
        AssertionError:
            If the subset name is not one of "train", "val", or "test".
            If the initial seed is not between 1 and 9999.
            If the sum of initial seed and count is greater than 10000.
    """

    SUBSET_OFFSETS = {"train": 10000, "val": 20000, "test": 30000}
    assert subset in SUBSET_OFFSETS, f"Unknown subset: {subset}, expected one of {list(SUBSET_OFFSETS.keys())}"
    assert 0 < init_seed < 10000, "Initial seed must be between 1 and 9999"
    assert init_seed + count < 10000, "Sum of initial seed and count must be less than 10000"

    return [SUBSET_OFFSETS[subset] + init_seed + i for i in range(count)]
