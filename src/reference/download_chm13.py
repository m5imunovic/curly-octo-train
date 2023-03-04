"""Extracts the chm13 reference genome into a fasta file per chromosome.

For simplicity, we ignore the last chromosome in the file.
"""

import multiprocessing as mp
import shlex
import shutil
import subprocess
from itertools import pairwise, repeat
from pathlib import Path

import requests
from omegaconf import OmegaConf
from torch_geometric.data import download_url, extract_gz

from reference.genome_generator import save_chr_to_fasta

# URLs of the chm13 reference genomes
chm13_urls = {
    "v1.1": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz",
    "v2.0": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
}


def download_if_not_exist(dowload_dir: Path, filename: str, url: str):
    download_path = dowload_dir / filename
    if not download_path.exists():
        print(f"Downloading reference from {url}...")
        requests.get(url, stream=True)


def download_chm13(dowload_dir: Path, url: str):
    filename = url.split("/")[-1]
    download_if_not_exist(dowload_dir, filename, url)


def grep_chromosome_offsets(chm13_path: Path) -> dict:
    # grep -n -e "^>" chm13.draft_v1.1.fasta
    cmd = shlex.join(["grep", "-n", "-o", "-E", "^>chr([0-9]*|[A-Z])", str(chm13_path)])
    print(cmd)
    grep_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    print(grep_output)
    offsets = {}
    for entry in grep_output.split():
        line_number, chromosome = entry.split(":>")
        offsets[chromosome] = int(line_number)

    return offsets


def get_chr_ranges(offsets: dict) -> dict:
    ranges = [(start, stop) for start, stop in pairwise(offsets.values())]
    # We don't have a range for the last chromosome
    chr_ranges = {chr: ranges[idx] for idx, chr in enumerate(offsets.keys()) if idx < len(ranges)}

    return chr_ranges


def generate_chr_fasta_files(
    chm13_path: Path, output_dir: Path, chr_ranges_entry: tuple
):  # , chr: str, chr_range: tuple):
    chromosome, (start, stop) = chr_ranges_entry
    cmd = shlex.join(["sed", "-n", f"{start},{stop-1}p", str(chm13_path)])
    print(cmd)
    sed_output = subprocess.check_output(cmd, shell=True, encoding="utf-8")
    chromosomes_dir = output_dir / "chromosomes"
    lines = sed_output.split('\n')
    chr_seq = "".join(lines[1:]).upper()
    save_chr_to_fasta(output_path=chromosomes_dir, chr_name=chromosome, chr_seq=chr_seq, multiline=False)


def main(cfg: OmegaConf):
    url = chm13_urls[cfg.chm13_version]
    download_dir = Path(cfg.output_dir)
    download_dir.mkdir(exist_ok=True, parents=True)
    fasta_gz_file = download_url(url=url, folder=download_dir)
    fasta_file = Path(fasta_gz_file).with_suffix("")
    if not fasta_file.exists():
        print(f"Extracting archive {fasta_gz_file} file...")
    extract_gz(path=fasta_gz_file, folder=download_dir)
    offsets = grep_chromosome_offsets(fasta_file)
    chr_ranges = get_chr_ranges(offsets)

    with mp.Pool(processes=cfg.threads) as pool:
        pool.starmap(generate_chr_fasta_files, zip(repeat(fasta_file), repeat(download_dir), chr_ranges.items()))

    fasta_file.unlink()
    if not cfg.debug:
        fasta_gz_file.unlink()

    shutil.move(str(download_dir), str(cfg.output_dir))


if __name__ == "__main__":
    # TODO: Use hydra.main and file config
    cfg = OmegaConf.create(
        {
            "chm13_version": "v2.0",
            "output_dir": "/home/${oc.env:USER}/ws/genomic/data/references/chm13_v2.0",
            "debug": True,
            "threads": 15,
        }
    )
    OmegaConf.resolve(cfg)
    main(cfg)
