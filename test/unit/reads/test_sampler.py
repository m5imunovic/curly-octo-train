import logging
from pathlib import Path

from omegaconf import DictConfig

from reads.rsimulator import READ_FILE
from reads.sampler import ReadSampler
from reads.simulate_reads import simulator_factory


def test_read_sampler_produces_expected_outputs(test_species_reads_root, tmpdir, caplog):
    cfg = DictConfig({"name": "sampler"})
    reads_simulator = simulator_factory("sampler", cfg)
    assert isinstance(reads_simulator, ReadSampler), "Constructed wrong simulator"

    test_species_reads = list(test_species_reads_root.glob("**/*.fastq"))
    output_path = Path(tmpdir)
    exec_args = {
        "reads": test_species_reads,
        "output_path": output_path,
        "probability": 0.25,
        "seed": 10,
    }

    with caplog.at_level(logging.DEBUG, logger="reads.sampler"):
        reads_simulator.run(**exec_args)

    assert (output_path / READ_FILE).exists()
