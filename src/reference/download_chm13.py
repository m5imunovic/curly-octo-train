"""Extracts the chm13 reference genome into a fasta file per chromosome. The results are stored in
output_dir/chromosomes directory.

For simplicity, we ignore the last chromosome (mitochondrial DNA) in the file.
"""

import logging
import multiprocessing as mp
import shlex
import subprocess
from itertools import pairwise, repeat
from pathlib import Path

import hydra
from omegaconf import OmegaConf
from torch_geometric.data import download_url, extract_gz

import utils.path_helpers as ph
import wandb  # TODO: resolve black and isort conflict on this line
from reference.genome_generator import save_chr_to_fasta

logger = logging.getLogger(__name__)
# URLs of the chm13 reference genomes
chm13_urls = {
    "v1.1": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz",
    "v2.0": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
}


def check_config(cfg: OmegaConf):
    """Checks the configuration for the download_chm13 script.

    Args:
        cfg (OmegaConf): app config
    """
    if cfg.chm13_version not in chm13_urls.keys():
        raise ValueError(f"Invalid chm13 version {cfg.chm13_version}. Valid values are: {chm13_urls.keys()}")
    if cfg.output_dir is None:
        raise ValueError("Output directory <output_dir> must be specified")
    # simple checks to capture config renames
    assert "threads" in cfg, "Missing threads parameter in config"
    assert "keep_metadata" in cfg, "Missing keep_metadata parameter in config"
    assert "store_artifacts" in cfg, "Missing store_artifacts parameter in config"


def grep_chromosome_offsets(chm13_path: Path) -> dict:
    """Find line number offsets for each chromosome in the chm13 fasta file. Grep is used to print the line number (-n)
    and the matched string (-o) for each chromosome. Extended grep (-E) is used to avoid backslash escaping of the
    regex.

    Args:
        chm13_path (Path): input file path

    Returns:
        dict: chromosome to offset mapping
    """
    cmd = shlex.join(["grep", "-n", "-o", "-E", "^>chr([0-9]*|[A-Z])", str(chm13_path)])
    logger.info(f"Running command: {cmd}")
    grep_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    logger.info(f"grep output: {grep_output}")
    offsets = {}
    for entry in grep_output.split():
        line_number, chromosome = entry.split(":>")
        offsets[chromosome] = int(line_number)

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


def generate_chr_fasta_files(chm13_path: Path, output_dir: Path, chr_ranges_entry: tuple):
    """Generate a fasta file for a single chromosome.

    Args:
        chm13_path (Path): input fasta file containing all chromosomes
        output_dir (Path): output location for the chromosome fasta files
        chr_ranges_entry (tuple): tuple of chromosome name and (start, stop) line numbers
    """
    chromosome, (start, stop) = chr_ranges_entry
    file_path = output_dir / f"{chromosome}.fasta"
    if file_path.exists():
        logger.info(f"The file path {str(file_path)} already exists, skipping...")
        return
    cmd = shlex.join(["sed", "-n", f"{start},{stop-1}p", str(chm13_path)])
    logger.info(f"Running command: {cmd}")
    sed_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    lines = sed_output.split("\n")
    chr_seq = "".join(lines[1:]).upper()
    save_chr_to_fasta(output_path=output_dir, chr_name=chromosome, chr_seq=chr_seq, multiline=False)


@hydra.main(version_base=None, config_path=str(ph.get_config_root() / "reference"), config_name="download_chm13")
def main(cfg: OmegaConf):
    url = chm13_urls[cfg.chm13_version]
    download_dir = Path(cfg.output_dir)
    download_dir.mkdir(exist_ok=True, parents=True)
    fasta_gz_file = download_url(url=url, folder=download_dir)
    fasta_file = Path(fasta_gz_file).with_suffix("")
    if not fasta_file.exists():
        logger.info(f"Extracting archive {fasta_gz_file} file...")
    extract_gz(path=fasta_gz_file, folder=download_dir)
    offsets = grep_chromosome_offsets(fasta_file)
    chr_ranges = get_chr_ranges(offsets)

    chromosomes_dir = download_dir / "chromosomes"
    chromosomes_dir.mkdir(exist_ok=True)
    with mp.Pool(processes=cfg.threads) as pool:
        pool.starmap(generate_chr_fasta_files, zip(repeat(fasta_file), repeat(chromosomes_dir), chr_ranges.items()))

    if cfg.store_artifacts:
        run = wandb.init(project="chm13", job_type="add-dataset")
        artifact = wandb.Artifact(name="chromosomes", type="dataset")
        artifact.add_dir(local_path=str(chromosomes_dir))
        run.log_artifact(artifact)

    fasta_file.unlink()
    if not cfg.keep_metadata:
        fasta_gz_file.unlink()


if __name__ == "__main__":
    main()
