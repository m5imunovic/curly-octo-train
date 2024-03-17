import logging
import shutil
import tempfile
from pathlib import Path

from typeguard import typechecked

from reads.rsimulator import READ_FILE, RSimulator
from reads.sample_reads import sample_reads_with_probability

logger = logging.Logger(__name__)
logger.setLevel(logging.INFO)


class ReadSampler(RSimulator):
    def _install(self, vendor_dir: Path):
        return vendor_dir

    def _sample_reads(self, reads: list[Path], output_path: Path, probability: float, seed: int) -> int:
        with tempfile.TemporaryDirectory() as staging_dir:
            sampled_file = Path(staging_dir) / READ_FILE
            total_selected = 0
            for read in reads:
                selected = sample_reads_with_probability(read, sampled_file, probability, seed)
                logger.info(f"Selected {selected} reads in {sampled_file} files.")
                total_selected += selected

            if not output_path.exists():
                output_path.mkdir(parents=True)
            shutil.copytree(staging_dir, output_path, dirs_exist_ok=True, copy_function=shutil.move)

            return total_selected

    @typechecked
    def run(self, reads: list[Path], output_path: Path, probability: float, seed: int, *args, **kwargs) -> bool:
        logger.info("Sampling reads...")
        selected = self._sample_reads(reads, output_path, probability, seed)
        logging.info(f"Selected total {selected} reads in {len(reads)} files.")
        return True
