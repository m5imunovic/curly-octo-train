import tempfile
from pathlib import Path
from unittest import mock

import pytest
from omegaconf import OmegaConf

from reads import reads_simulator
from reads.reads_simulator import PbSim3
from reads.reads_simulator import run as run_reads_simulator

SEED = 1


def reads_simulator_cfg(overwrite: bool, tmp_dir: Path):
    return OmegaConf.create(
        {
            "reads": {
                "name": "pbsim3",
                "overwrite": overwrite,
                "params": {"long": {"depth": 30, "seed": SEED}},
                "request": {"chr1": 2},
            },
            "paths": {
                "ref_dir": str(tmp_dir / "ref"),
                "reads_dir": str(tmp_dir / "simulated"),
            },
            "date_mm_dd": "01_01",
            "species_name": {"name": "test_species"},
            "experiment": "test_experiment",
            "seed": SEED,
        }
    )


@pytest.mark.skip(reason="Not compatible with profile simulation")
@mock.patch("random.randint", return_value=22)
def test_pbmsim3_construct_exe_cmd(mock_randint, test_reads_root):
    cfg1 = OmegaConf.create(
        {
            "name": "pbsim3",
            "params": {
                "long": {
                    "depth": 30,
                    "seed": SEED,
                },
                "pattern": None,
            },
        }
    )

    cfg2 = OmegaConf.create({"name": "pbsim3", "params": {}})

    with tempfile.TemporaryDirectory() as tmp_dir:
        vendor_dir = Path(tmp_dir) / "vendor"
        vendor_dir.mkdir(exist_ok=True)
        simulator_root = vendor_dir / "pbsim3"
        simulator_root.mkdir(exist_ok=True)

        dummy_ref = test_reads_root / "dummy" / "dummy.fa"
        save_dir = Path(tmp_dir) / "save"
        prefix = "1"
        pbsim = reads_simulator.PbSim3(cfg1, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)
        seed = mock_randint.return_value + int(prefix)

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")} --depth 30 --prefix {prefix} {str(dummy_ref)} --seed {seed}',
            f"rm {prefix}_0001.ref",
        ]

        assert len(cmd) == 2
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]

        pbsim = reads_simulator.PbSim3(cfg2, vendor_dir)
        cmd = pbsim._construct_exec_cmd(dummy_ref, save_dir, prefix)
        assert len(cmd) == 2

        expected_cmd = [
            f'{str(simulator_root / "src/pbsim")}  --prefix {prefix} {str(dummy_ref)} --seed {seed}',
            f"rm {prefix}_0001.ref",
        ]
        assert cmd[0] == expected_cmd[0]
        assert cmd[1] == expected_cmd[1]


@mock.patch.object(PbSim3, "run")
def test_reads_simulator_overwrite_data(mock_run, tmp_path):
    with tempfile.TemporaryDirectory() as tmp_dir:
        cfg = reads_simulator_cfg(overwrite=True, tmp_dir=Path(tmp_dir))
        output_path = Path(cfg.paths.reads_dir) / cfg.species_name["name"] / cfg.date_mm_dd / f"S{cfg.seed}"
        kwargs = {
            "ref_root": tmp_path / cfg.species_name["name"],
            "chr_request": dict(cfg.reads.request),
            "simulated_species_path": output_path,
        }
        mock_run.return_value = True
        output_path.mkdir(parents=True, exist_ok=True)
        run_reads_simulator(cfg, **kwargs)
        assert mock_run.call_count == 1

        cfg = reads_simulator_cfg(overwrite=False, tmp_dir=Path(tmp_dir))
        with pytest.raises(FileExistsError):
            run_reads_simulator(cfg, **kwargs)
