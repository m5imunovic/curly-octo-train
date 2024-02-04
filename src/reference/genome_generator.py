import json
import logging
import shutil
from pathlib import Path

import reference.reference_utils as ru
from reference.chm13 import get_chm13_reference
from reference.random_genome import get_random_reference

logger = logging.Logger(__name__)


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
        species_path.mkdir(parents=True)

        if ru.is_random_reference(species):
            genome = get_random_reference(dict(species.chromosomes), species.seed, species.gc_content)
            ru.save_genome_to_fasta(ru.ref_chromosomes_path(species_path), genome, multiline=False)
            msg = f"Generated random reference genome for species `{species.name}` at:  \n{reference_root}"
            logger.info(msg)
        else:
            # TODO: save url from which reference was downloaded to species_info.json
            get_chm13_reference(species_path, species)

        with open(species_path / "species_info.json", "w") as handle:
            species_info = ru.create_species_info(species)
            json.dump(species_info, handle, indent=4)

    except Exception as e:
        shutil.rmtree(species_path)
        msg = f"Error while generating reference genome for species `{species.name}`: {e}"
        return None

    return species_path
