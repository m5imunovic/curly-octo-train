import json
import logging
import shutil
import subprocess
from copy import deepcopy
from pathlib import Path

from omegaconf import OmegaConf
from typeguard import typechecked

from asm import assembler
from utils.io_utils import compose_cmd_params

logger = logging.getLogger(__name__)


class LaJolla(assembler.Assembler):
    @typechecked
    def _install(self, vendor_dir: Path) -> Path:
        assembler_root = (vendor_dir / "LJA").absolute()
        if not assembler_root.exists():
            try:
                logger.info("SETUP::generate:: Download La Jolla Assembler")
                subprocess.run(
                    "git clone https://github.com/AntonBankevich/LJA.git", shell=True, cwd=str(vendor_dir.absolute())
                )
                subprocess.run("git fetch", shell=True, cwd=str(assembler_root))
                subprocess.run("git checkout -t origin/anton_development", shell=True, cwd=str(assembler_root))
                subprocess.run("cmake .", shell=True, cwd=str(assembler_root))
                subprocess.run("make -j 8 lja", shell=True, cwd=str(assembler_root))
                subprocess.run("make -j 8 jumboDBG", shell=True, cwd=str(assembler_root))
                subprocess.run("make -j 8 align_and_print", shell=True, cwd=str(assembler_root))
            except Exception as ex:
                logger.info(f"SETUP::generate:: Error: {ex}")
                shutil.rmtree(assembler_root)
                raise ex

        return assembler_root / "bin"

    @typechecked
    def _construct_eval_cmds(self, output_path: Path, reads_cmd_params: str) -> Path | None:
        cmds = {}
        if self.cfg["full_asm"]:
            full_asm_subdir = self.cfg["full_asm"]
            params = OmegaConf.to_container(self.cfg["params"])
            params["long"].pop("reference", None)
            params["append"] = params["append"].replace("--compress", "")
            params["short"].pop("K", None)
            k = params["short"].pop("k", None)

            eval_cmds_dir = output_path / full_asm_subdir
            eval_cmds_dir.mkdir(exist_ok=True, parents=True)
            eval_cmds_file = eval_cmds_dir / "full_asm.json"
            read_files = reads_cmd_params.split()
            for idx, read_file in enumerate(read_files):
                if idx % 2:
                    rf = Path(read_file)
                    rf = Path("${eval_dir}/") / rf.relative_to(rf.parent.parent.parent)
                    read_files[idx] = str(rf)
            reads_cmd_params = " ".join(read_files)

            output_path = Path("${eval_dir}/") / output_path.relative_to(output_path.parent.parent)
            full_asm_path = output_path / full_asm_subdir
            params["short"]["o"] = str(full_asm_path)
            params["append"] = reads_cmd_params
            cmds["lja"] = {"params": deepcopy(params)}

            params = {
                "short": {
                    "k": k,
                },
                "long": {
                    "dbg": str(full_asm_path / "00_CoverageBasedCorrection" / "initial_dbg.gfa"),
                    "paths": str(full_asm_path / "00_CoverageBasedCorrection" / "corrected_reads.fasta"),
                },
            }

            params["short"]["o"] = str(output_path / "eval_00")

            cmds["align_and_print"] = {}
            # command for generating alignments
            cmds["align_and_print"]["eval_00"] = deepcopy(params)

            params["short"]["o"] = str(output_path / "eval_01")
            params["long"]["paths"] = str(full_asm_path / "01_TopologyBasedCorrection" / "corrected_reads.fasta")
            # command for generating alignments
            cmds["align_and_print"]["eval_01"] = deepcopy(params)

            with open(eval_cmds_file, "w") as f:
                json.dump(cmds, f, indent=4)
            return eval_cmds_file

        return None

    @typechecked
    def _construct_exec_cmd(
        self, genome: list[Path], reads: list[Path], output_path: Path
    ) -> tuple[list[str], Path | None]:
        assert "params" in self.cfg, "params must be specified in config"

        if genome:
            self.cfg.params.long.reference = [str(p) for p in genome]
        self.cfg.params.short.o = str(output_path)
        option_params = compose_cmd_params(self.cfg["params"])

        exec = "lja" if "exec" not in self.cfg else self.cfg["exec"]
        asm_executable = self.assembler_root / exec
        # TODO: handle reads through config
        reads_cmd_params = " ".join([f"--reads {str(read_file)}" for read_file in reads])
        cmds = [f"{asm_executable} {option_params} {reads_cmd_params}"]

        eval_cmds_path = self._construct_eval_cmds(output_path=output_path, reads_cmd_params=reads_cmd_params)

        return cmds, eval_cmds_path

    @typechecked
    def run(self, genome: list[Path], reads: list[Path], output_path: Path, *args, **kwargs):
        commands, eval_cmds_path = self._construct_exec_cmd(genome, reads, output_path)
        logger.info(f"Save eval commands to {eval_cmds_path}")

        for cmd in commands:
            logger.info(f"RUN::assembler::\n{cmd}")
            subprocess.run(cmd, shell=True, cwd=output_path)

    @typechecked
    def pre_assembly_step(self, output_path: Path, *args, **kwargs):
        if output_path.exists():
            logger.info("PRE::assembler:: Warning, existing assembly files, overwriting them")
            shutil.rmtree(output_path)
        output_path.mkdir(parents=True)

    def post_assembly_step(self, output_path: Path, *args, **kwargs):
        # Clean up the files
        if "keep" in self.cfg and self.cfg["keep"] is not None:
            keep_files = OmegaConf.to_container(self.cfg["keep"])
            for path in output_path.glob("**/*"):
                if not path.is_dir():
                    if self.cfg["full_asm"]:
                        if "full_asm" in str(path):
                            continue
                    if path.name not in keep_files:
                        path.unlink()
