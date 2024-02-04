from pathlib import Path

from experiment.schedule import create_sequencing_jobs


def test_create_sequencing_jobs(tmpdir, test_scenario, test_data_reference):
    reference_paths = {"test_species.yaml": test_data_reference / "test_species/S_0001"}
    jobs = create_sequencing_jobs(test_scenario, sequencing_root=Path(tmpdir), reference_paths=reference_paths)
    assert len(jobs) == 3
    print(jobs)

    job = jobs[0]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 10001
    assert job["simulated_reads_path"].name == "test_species_S10001_chr1"

    job = jobs[1]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["seed"] == 10002
    assert job["simulated_reads_path"].name == "test_species_S10002_chr1"

    job = jobs[2]
    assert job["genome"][0].name == "chr1.fasta"
    assert job["genome"][1].name == "chr2.fasta"
    assert job["seed"] == 10001
    assert job["simulated_reads_path"].name == "test_species_S10001_chr1_chr2"
