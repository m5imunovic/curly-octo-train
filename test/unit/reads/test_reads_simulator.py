import filecmp
import tempfile
from pathlib import Path

import pytest

from reads import reads_simulator
from reads.reads_simulator import simulator_factory
from reads.simulate_reads import SampleProfile

SEED = 3


@pytest.fixture(scope="module")
def fake_vendor_root():
    # TODO: use pyfakefs for this
    with tempfile.TemporaryDirectory() as tmp_dir:
        vendor_dir = Path(tmp_dir) / "vendor"
        vendor_dir.mkdir(exist_ok=True)
        simulator_root = vendor_dir / "pbsim3"
        simulator_root.mkdir(exist_ok=True)
        # PbSim3 will expects vendor_dir as input argument and that pbsim3 folder exists in there
        yield vendor_dir


@pytest.fixture(scope="module")
def fake_sample_profile():
    profile = SampleProfile(
        genome=[Path("/fakereferences/random_species/S_0100/chromosomes/chr1.fasta")],
        seed=SEED,
    )
    return profile


@pytest.fixture(scope="function")
def test_species_genome(test_data_reference) -> list[Path]:
    chr_path = test_data_reference / "test_species" / "S_0001" / "chromosomes"
    return list(chr_path.glob("*.fasta"))


def test_generator_produces_expected_output(test_pbsim3_pf_cfg, test_species_genome, test_data_reads, tmpdir):
    simulator = simulator_factory("pbsim3", test_pbsim3_pf_cfg)
    simulator.run(genome=test_species_genome, simulated_reads_path=Path(tmpdir))

    simulated_fastq = tmpdir / "sim_0001.fastq"
    assert simulated_fastq.exists(), "Simulated fastq file does not exist"
    expected_fastq = test_data_reads / "simulated" / "sim_0001.fastq"
    assert filecmp.cmp(simulated_fastq, expected_fastq), "Simulated fastq file does not match expected file"


def test_simulator_generates_expected_commands(test_pbsim3_pf_cfg, fake_vendor_root, fake_sample_profile):
    profile_root = Path("/home/fake/pf/")
    cfg = test_pbsim3_pf_cfg
    # override parameters coming from the profile
    cfg.params.long.seed = fake_sample_profile.seed
    cfg.profile.path = str(profile_root)

    pbsim = reads_simulator.PbSim3(cfg, fake_vendor_root)
    cmd = pbsim.construct_exec_cmd(fake_sample_profile.genome)

    simulator_exe = fake_vendor_root / "pbsim3/src/pbsim"
    profile_path = pbsim.construct_sample_path(profile_root, cfg.params.long["sample-profile-id"], ".fastq")
    stats_path = pbsim.construct_sample_path(profile_root, cfg.params.long["sample-profile-id"], ".stats")
    expected_cmds = [
        f"ln {profile_path} {profile_path.name}",
        f"ln {stats_path} {stats_path.name}",
        f"{simulator_exe} --depth 25 --strategy wgs --method sample --sample-profile-id pf1 --seed {SEED} --prefix sim --genome {fake_sample_profile.genome[0]}",
        f"rm {cfg.params.long.prefix}_0001.maf",
        f"rm {cfg.params.long.prefix}_0001.ref",
        f"rm {profile_path.name}",
        f"rm {stats_path.name}",
    ]

    for c, e in zip(cmd, expected_cmds):
        assert c == e, f"{c} != {e}"
