import json
import logging
import shutil
from copy import deepcopy
from pathlib import Path

from torch_geometric.data import download_url

import reference.reference_utils as ru
from reference.bed.bed_parser import parse_bed_file, merge_bed_regions, pretty_print_regions

logger = logging.getLogger(__name__)


def run(cfg, neighborhood: int = 100000, gap: int = 20000):
    reference_root = Path(cfg.paths.ref_dir)
    species = cfg.reference.species
    species_path = ru.get_species_root(Path(reference_root), species)

    bed_species = deepcopy(species)
    bed_species["release"] = "bed_" + bed_species["release"]

    bed_species_path = ru.get_species_root(Path(reference_root), bed_species)
    bed_species_path.mkdir(parents=True, exist_ok=True)

    # TODO: check bed species exists

    try:
        # download bed file
        bed_file = download_url(url=cfg.reference.bed.url, folder=bed_species_path)
        bed_entries = merge_bed_regions(parse_bed_file(bed_file), neighborhood=neighborhood, gap=gap)

        chromosome_paths = ru.query_chr_paths(ru.ref_chromosomes_path(species_path))

        # filter out chromosome paths which don't have bed input
        missing_beds = set(chromosome_paths) - set(bed_entries)
        chromosome_paths = {
            chr_name: chr_path for chr_name, chr_path in chromosome_paths.items() if chr_name not in missing_beds
        }

        # filter out bed ranges for which we did not found chromosome and issue warning
        missing_chrs = set(bed_entries) - set(chromosome_paths)
        if len(missing_chrs) > 0:
            logger.warning(f"Missing chromosomes {','.join(missing_chrs)} referenced in bed file")
        chromosome_paths = {
            chr_name: chr_path for chr_name, chr_path in chromosome_paths.items() if chr_name not in missing_chrs
        }

        extracted_meta = {}
        for chromosome_name, chromosome_path in chromosome_paths.items():
            if bed_entries[chromosome_name] is None:
                logger.warning(f"Skipping {chromosome_name} as the range is not well defined")
                continue
            subreference, interval = ru.extract_chromosome_range(chromosome_path, bed_entries[chromosome_name])
            # save into bed reference
            ru.save_genome_to_fasta(ru.ref_chromosomes_path(bed_species_path), subreference, multiline=False)
            extracted_meta[chromosome_name] = interval

        pretty_print_regions(extracted_meta, bed_species_path / "regions.tsv")

        with open(bed_species_path / "species_info.json", "w") as handle:
            species_info = ru.create_species_info(bed_species)
            json.dump(species_info, handle, indent=4)

    except Exception as e:
        shutil.rmtree(bed_species_path)
        msg = f"Error while generating reference for `{bed_species.name}` from bed: {e}"
        logger.error(msg)
        return None

    return bed_species_path
