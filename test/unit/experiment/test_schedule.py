from pathlib import Path
from unittest import mock

from experiment.schedule import create_sequencing_jobs, run


def test_create_sequencing_jobs(tmpdir, test_scenario_sim, test_data_reference):
    reference_paths = {"test_species.yaml": test_data_reference / "test_species/S_0001"}
    staging_root = Path(tmpdir)
    jobs = create_sequencing_jobs(test_scenario_sim, staging_root=staging_root, reference_paths=reference_paths)
    assert len(jobs) == 3

    job = jobs[0]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 30001
    assert job["output_path"].relative_to(staging_root) == Path("0/reads")

    job = jobs[1]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 30002
    assert job["output_path"].relative_to(staging_root) == Path("1/reads")

    job = jobs[2]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["genome"][1].name == "chr2.fasta"
    assert job["seed"] == 30001
    assert job["output_path"].relative_to(staging_root) == Path("2/reads")


def test_schedule_run_sim_creates_expected_outputs(
    test_cfg_root, test_experiment_sim_cfg, test_scenarios_root, tmpdir
):
    test_experiment_sim_cfg.paths.exp_dir = str(tmpdir)
    test_experiment_sim_cfg.paths.datasets_dir = str(tmpdir)
    with mock.patch("experiment.experiment_utils.get_scenario_root", return_value=test_scenarios_root), mock.patch(
        "utils.path_helpers.get_config_root", return_value=test_cfg_root
    ):
        run(test_experiment_sim_cfg)

    assert (Path(tmpdir) / "unittest_dataset").exists()
    assert len(list((Path(tmpdir) / "unittest_dataset" / "test").glob("**/*.pt"))) == 3
    eval_path = Path(tmpdir) / "unittest_dataset" / "eval"
    assert eval_path.exists()
    for sample_idx in ["0", "1", "2"]:
        # do hard-coded checking, actually should be cross-checked with experiment.keep option
        assert eval_path / sample_idx / "reads" / "sim_0001.fastq"
        assert eval_path / sample_idx / "assemblies" / "mult.info"
        assert eval_path / sample_idx / "graph" / "debug" / "0.idmap"


def test_schedule_run_sample_creates_expected_outputs(
    test_cfg_root, test_experiment_sampler_cfg, test_scenarios_root, tmpdir
):
    test_experiment_sampler_cfg.paths.exp_dir = str(tmpdir)
    test_experiment_sampler_cfg.paths.datasets_dir = str(tmpdir)
    with mock.patch("experiment.experiment_utils.get_scenario_root", return_value=test_scenarios_root), mock.patch(
        "utils.path_helpers.get_config_root", return_value=test_cfg_root
    ):
        run(test_experiment_sampler_cfg)

    assert (Path(tmpdir) / "unittest_dataset").exists()

    assert len(list((Path(tmpdir) / "unittest_dataset" / "train" / "raw").iterdir())) == 3
