from pathlib import Path
from unittest import mock

from experiment.schedule import create_sequencing_jobs, run


def test_create_sequencing_jobs(tmpdir, test_scenario, test_data_reference):
    reference_paths = {"test_species.yaml": test_data_reference / "test_species/S_0001"}
    jobs = create_sequencing_jobs(test_scenario, sequencing_root=Path(tmpdir), reference_paths=reference_paths)
    assert len(jobs) == 3
    print(jobs)

    job = jobs[0]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 10001
    assert job["output_path"].name == "test_species_S10001_chr1"

    job = jobs[1]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 10002
    assert job["output_path"].name == "test_species_S10002_chr1"

    job = jobs[2]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["genome"][1].name == "chr2.fasta"
    assert job["seed"] == 10001
    assert job["output_path"].name == "test_species_S10001_chr1_chr2"


def test_schedule_run_creates_expected_outputs(test_cfg_root, test_experiment_cfg, test_scenarios_root, tmpdir):
    test_experiment_cfg.paths.exp_dir = str(tmpdir)
    with mock.patch("experiment.experiment_utils.get_scenario_root", return_value=test_scenarios_root), mock.patch(
        "utils.path_helpers.get_config_root", return_value=test_cfg_root
    ):
        run(test_experiment_cfg)

    assert (Path(tmpdir) / "unittest_dataset").exists()
