import json
import logging
import shutil
from pathlib import Path

from omegaconf import DictConfig, OmegaConf

import experiment.experiment_utils as eu
import reference.reference_utils as ru
import utils.path_helpers as ph
from reference.bed_genome_generator import run as bed_run
from reference.chm13 import get_chm13_reference
from reference.hg002 import get_hg002_reference

# from reference.pan_troglodytes import get_pan_tro3_reference
from reference.random_genome import get_random_reference
from reference.repeat_sequence import get_repeat_sequence

logger = logging.getLogger(__name__)


def run(cfg):
    reference_root = Path(cfg.paths.ref_dir)
    species = cfg.reference.species
    species_path = ru.get_species_root(Path(reference_root), species)

    # TODO: this function should accept species_path
    if ru.check_reference_exists(reference_root, species):
        msg = f"Reference genome for species `{species.name}` already exists at:  \n{reference_root}"
        logger.info(msg)
        return species_path

    try:
        species_path = ru.get_species_root(Path(reference_root), species)
        species_path.mkdir(parents=True, exist_ok=True)

        if ru.is_random_reference(species):
            genome = get_random_reference(dict(species.chromosomes), species.seed, species.gc_content)
            ru.save_genome_to_fasta(ru.ref_chromosomes_path(species_path), genome, multiline=False)
            msg = f"Generated random reference genome for species `{species.name}` at:  \n{reference_root}"
            logger.info(msg)
        elif ru.is_chm13_reference(species):
            # TODO: save url from which reference was downloaded to species_info.json
            get_chm13_reference(species_path, species)
        elif ru.is_hg002_reference(species):
            get_hg002_reference(species_path, species)
        # elif ru.is_mPanTro3_reference(species):
        #    get_pan_tro3_reference(species_path, species)
        elif ru.is_repeat_reference(species):
            get_repeat_sequence(cfg.reference.species, species_path)
        else:
            raise ValueError(f"Unsupported reference {species}")

        with open(species_path / "species_info.json", "w") as handle:
            species_info = ru.create_species_info(species)
            json.dump(species_info, handle, indent=4)

    except Exception as e:
        shutil.rmtree(species_path)
        msg = f"Error while generating reference genome for species `{species.name}`: {e}"
        logger.error(msg)
        return None

    return species_path


def ensure_references_exist(cfg: DictConfig, species_defs: set[str]) -> dict | None:
    """Ensures that all species references exist for the given scenario.

    Args:
        species_defs (set[str]): The set of species definitions described in respective config file names.

    Returns:
        bool: True if all species references exist, False otherwise.
    """
    config_root = ph.get_config_root()
    try:
        eu.ensure_species_def_exists(config_root, species_defs)
    except AssertionError as e:
        logger.error(f"Species definition files error:\n {e}")
        return False

    # TODO: this is a candidate for parallelization
    # TODO: this is a bit ugly now how it is handled, we should use hydra compose to load the files
    species_paths = {}
    species_config_root = config_root / "reference" / "species"
    for species_def in species_defs:
        cfg_species = OmegaConf.load(species_config_root / f"{species_def}")
        defaults = cfg_species.pop("defaults", None)
        if defaults is not None:
            # we assume only one defaults entry and that it exists
            cfg_defaults = OmegaConf.load(species_config_root / f"{defaults[0]}.yaml")
            cfg_species = OmegaConf.merge(cfg_defaults, cfg_species)

        cfg.reference.species = cfg_species
        reference_path = run(cfg)
        if not reference_path:
            reference_path = bed_run(cfg)
            if not reference_path:
                return None
        species_paths[species_def] = reference_path

    return species_paths
