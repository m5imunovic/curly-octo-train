"""Extracts the HG002 reference genome into a fasta file per chromosome."""

import logging
import multiprocessing as mp
import shlex
import subprocess
from itertools import pairwise, repeat
from pathlib import Path

import wandb
from omegaconf import DictConfig
from torch_geometric.data import download_url, extract_gz

import reference.reference_utils as ru

logger = logging.getLogger(__name__)
# URLs of the HG002 reference genomes
hg002_urls = {
    "hg002_v1_0_1": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz"
}


def check_config(cfg: DictConfig):
    """Checks the configuration for the hg002.

    Args:
        cfg (DictConfig): species config
    """
    if cfg.release not in hg002_urls.keys():
        raise ValueError(f"Invalid hg002 version {cfg.release}. Valid values are: {hg002_urls.keys()}")


def grep_chromosome_offsets(ref_fasta_path: Path) -> dict:
    """Find line number offsets for each chromosome in the reference fasta file. Grep is used to print the line number
    (-n) and the matched string (-o) for each chromosome. Extended grep (-E) is used to avoid backslash escaping of the
    regex.

    Args:
        ref_fasta_path (Path): input file path

    Returns:
        dict: chromosome to offset mapping
    """
    cmd = shlex.join(["grep", "-n", "-o", "-E", "^>chr([0-9]|_|[A-Z])*", str(ref_fasta_path)])
    logger.info(f"Running command: {cmd}")
    grep_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    logger.info(f"grep output: {grep_output}")
    offsets = {}
    for entry in grep_output.split():
        line_number, chromosome = entry.split(":>")
        offsets[chromosome] = int(line_number)

    cmd = shlex.join(["wc", "-l", str(ref_fasta_path)])
    logger.info(f"Running command: {cmd}")
    wc_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    line_count = int(wc_output.split()[0])
    # artificial offset to denote the end for last chromosome in file
    offsets["END"] = line_count + 1

    return offsets


def get_chr_ranges(offsets: dict) -> dict:
    """Get the start and stop line numbers for each chromosome.

    Args:
        offsets (dict): chromosome to offset mapping

    Returns:
        dict: chromosome to (start, stop) mapping
    """
    ranges = [(start, stop) for start, stop in pairwise(offsets.values())]
    # We don't have a range for the last chromosome as end info is missing in offsets
    chr_ranges = {chr: ranges[idx] for idx, chr in enumerate(offsets.keys()) if idx < len(ranges)}

    return chr_ranges


def generate_chr_fasta_files(ref_fasta_path: Path, output_dir: Path, chr_ranges_entry: tuple):
    """Generate a fasta file for a single chromosome.

    Args:
        ref_fasta_path (Path): input fasta file containing all chromosomes
        output_dir (Path): output location for the chromosome fasta files
        chr_ranges_entry (tuple): tuple of chromosome name and (start, stop) line numbers
    """
    chromosome, (start, stop) = chr_ranges_entry
    file_path = output_dir / f"{chromosome}.fasta"
    if file_path.exists():
        logger.info(f"The file path {str(file_path)} already exists, skipping...")
        return
    cmd = shlex.join(["sed", "-n", f"{start},{stop-1}p", str(ref_fasta_path)])
    logger.info(f"Running command: {cmd}")
    sed_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    lines = sed_output.split("\n")
    chr_seq = "".join(lines[1:]).upper()
    ru.save_chr_to_fasta(output_path=output_dir, chr_name=chromosome, chr_seq=chr_seq, multiline=False)


def get_hg002_reference(species_root: Path, species: DictConfig, log: bool = True):
    """Download HG002 reference archive and unpack it in appropriate form."""
    url = hg002_urls[species.release]
    download_dir = species_root
    fasta_gz_file = download_url(url=url, folder=download_dir, log=log)
    fasta_file = Path(fasta_gz_file).with_suffix("")
    if not fasta_file.exists():
        logger.info(f"Extracting archive {fasta_gz_file} file...")
        extract_gz(path=fasta_gz_file, folder=download_dir, log=log)
        print("Extracted archive...")

    offsets = grep_chromosome_offsets(fasta_file)
    chr_ranges = get_chr_ranges(offsets)
    chromosomes_dir = ru.ref_chromosomes_path(download_dir)
    chromosomes_dir.mkdir(exist_ok=True)
    threads = min(mp.cpu_count() - 1, len(chr_ranges))
    with mp.Pool(processes=threads) as pool:
        pool.starmap(generate_chr_fasta_files, zip(repeat(fasta_file), repeat(chromosomes_dir), chr_ranges.items()))

    # TODO: move store_artifacts to reference config
    if species.store_artifacts:
        run = wandb.init(project="hg002", job_type="add-dataset")
        artifact = wandb.Artifact(name="chromosomes", type="dataset")
        artifact.add_dir(local_path=str(chromosomes_dir))
        run.log_artifact(artifact)

    fasta_file.unlink()
    # TODO move keep_metadata to reference config
    if not species.keep_metadata:
        fasta_gz_file.unlink()

    return url
