"""Evaluates the corrected reads produced by La Jolla Assembler.

More specifically, it evaluates the corrected reads produced by the second step of the pipeline (i.e. the topology-
based correction).
"""

import json
import logging
import os
import subprocess
from pathlib import Path
from typing import MutableSequence

import hydra
from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from asm.collect_evaluation import get_eval_table
from asm.mult_info_parser import parse_mult_info, partition_mult_info_edges
from graph.tools.gfa_parser import parse_gfa
from utils.io_utils import compose_cmd_params

logger = logging.getLogger(__name__)


def parse_alignments_entry(read_id: str, edge_ids: str):
    """First line is read ID, second line is list of edge IDs.

    Returns a tuple of (read_id, edge_ids)
    """
    read_id = read_id.strip()[1:]  # remove starting '>'
    edge_ids = edge_ids.strip().split()

    return read_id, edge_ids


@typechecked
def parse_alignments(alignments_path: Path) -> set:
    with open(alignments_path) as f:
        correct_reads = set()
        for read_line in f:
            _, edge_ids = parse_alignments_entry(read_line, next(f))
            correct_reads.update(edge_ids)

    return correct_reads


@typechecked
def get_correct_edges(alignments_path: Path) -> set:
    correct_edges_with_rc = parse_alignments(alignments_path)
    return correct_edges_with_rc


def map_lja_ids_to_jumbo_ids(initial_dbg_path, graph_dbg_path):
    logger.info("Get the LJA to JumboDBG mapping...")
    segments, _ = parse_gfa(initial_dbg_path)
    with open(graph_dbg_path / "0.hashmap") as f:
        jumbo_hashmap = json.load(f)

    lja_jumbo_map = {}
    count_misses = 0
    logger.info("Parsed LJA graph, perform mapping...")
    for seg_id, metadata in segments.items():
        hash = metadata["hash"]
        if hash in jumbo_hashmap:
            lja_jumbo_map[seg_id] = jumbo_hashmap[hash]
        else:
            count_misses += 1

    logger.info(f"Could not find {count_misses} LJA hashes in Jumbo hash map")

    return lja_jumbo_map


@typechecked
def get_confusion_matrix(mult_info_path: Path, alignments_path: Path, mapped_ids) -> tuple:
    """Requires the following files:

    - alignments.txt
    - mult.info
    """

    # get ground truth
    mult_info = parse_mult_info(mult_info_path)
    correct_edges_gt, incorrect_edges_gt = partition_mult_info_edges(mult_info)
    logger.info("Ground truth:")
    logger.info(f"\tCorrect edges: {len(correct_edges_gt)}")
    logger.info(f"\tIncorrect edges: {len(incorrect_edges_gt)}")

    # get mowerDBG edge assignments
    correct_edges = get_correct_edges(alignments_path)
    correct_edges = {mapped_ids[sid] for sid in correct_edges if sid in mapped_ids}
    incorrect_edges = (correct_edges_gt | incorrect_edges_gt) - correct_edges
    logger.info("MowerDBG classified:")
    logger.info(f"\t{len(correct_edges)} as correct edges.")
    logger.info(f"\t{len(incorrect_edges)} as incorrect edges.")

    # true positives
    tp = len(correct_edges_gt & correct_edges)
    tn = len(incorrect_edges_gt & incorrect_edges)
    fp = len(incorrect_edges_gt & correct_edges)
    fn = len(correct_edges_gt & incorrect_edges)

    return tp, tn, fp, fn


@typechecked
def calculate_metrics(tp: int, tn: int, fp: int, fn: int) -> dict[str, float]:
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    f1 = 2 * precision * recall / (precision + recall)
    return {
        "Precision": precision,
        "Recall": recall,
        "Accuracy": accuracy,
        "F1": f1,
    }


def evaluate_la_jolla(mult_info_path: Path, alignments_path: Path, mapped_ids: dict) -> dict:
    tp, tn, fp, fn = get_confusion_matrix(mult_info_path, alignments_path, mapped_ids)
    metrics = calculate_metrics(tp, tn, fp, fn)

    evaluations = {
        "Mult info path": str(mult_info_path),
        "Alignments path": str(alignments_path),
        "Confusion matrix": {
            "tp fn": [tp, fn],
            "fp tn": [fp, tn],
        },
        "metrics": metrics,
    }

    return evaluations


class PathDecoder(json.JSONDecoder):
    def __init__(self, eval_dir: Path):
        super().__init__()
        self.eval_dir = eval_dir

    def decode(self, s):
        decoded = super().decode(s)
        return self.replace_placeholder(decoded)

    def replace_placeholder(self, data):
        if isinstance(data, dict):
            for key, value in data.items():
                if isinstance(value, str):
                    data[key] = value.replace("${eval_dir}", self.eval_dir)
                elif isinstance(value, (dict, list)):
                    self.replace_placeholder(value)
        elif isinstance(data, list):
            for index, item in enumerate(data):
                if isinstance(item, str):
                    data[index] = item.replace("${eval_dir}", self.eval_dir)
                elif isinstance(item, (dict, list)):
                    self.replace_placeholder(item)
        return data


@typechecked
def construct_eval_commands(
    lja_bin_path: Path,
    eval_cmds_path: Path,
    skip_cmds: MutableSequence,
    eval_stages: MutableSequence,
    decoder_path: str,
) -> list[str]:
    eval_cmds = {}
    with open(eval_cmds_path) as f:
        json_data = f.read()
        decoder = PathDecoder(decoder_path)
        eval_cmds = decoder.decode(json_data)
    if not eval_cmds:
        raise ValueError("No commands found in eval_cmds.json file")

    cmds = []
    for lja_cmd in eval_cmds:
        if lja_cmd not in skip_cmds:
            if lja_cmd == "lja":
                params = eval_cmds[lja_cmd]["params"]
                cmds.append(f"{str(lja_bin_path / lja_cmd)} {compose_cmd_params(params)}")
            elif lja_cmd == "align_and_print":
                for stage in eval_cmds[lja_cmd]:
                    if stage in eval_stages:
                        params = eval_cmds[lja_cmd][stage]
                        cmds.append(f"{str(lja_bin_path / lja_cmd)} {compose_cmd_params(params)}")

            else:
                raise ValueError(f"Unknown command {lja_cmd}")

    return cmds


def eval_lja(cfg: DictConfig, subdir: str):
    output_path_eval = Path(cfg.eval_path) / subdir / "eval_results"
    if output_path_eval.exists():
        logger.warning(f"{output_path_eval} already exists, skipping evaluation for {subdir=}")
        return

    # First check that required directories and files are available
    asm_path = Path(cfg.eval_path) / subdir / "assemblies"
    mult_info_path = Path(asm_path) / "mult.info"
    assert mult_info_path.exists(), f"Mult info file {mult_info_path} does not exist"
    eval_cmds_path = asm_path / cfg.full_asm_subdir / cfg.eval_cmds_path
    assert eval_cmds_path.exists(), f"Evaluation commands path {eval_cmds_path} does not exist"
    graph_dbg_path = Path(cfg.eval_path) / subdir / "graph" / "debug"
    assert (graph_dbg_path / "0.idmap").exists()
    assert (graph_dbg_path / "0.hashmap").exists()

    logger.info("Constructing execution commands...")
    cmds = construct_eval_commands(
        lja_bin_path=Path(cfg.lja_bin_path),
        eval_cmds_path=eval_cmds_path,
        skip_cmds=cfg.skip_cmds,
        eval_stages=cfg.eval_stages,
        decoder_path=cfg.eval_path,
    )

    for cmd in cmds:
        logger.info(f"Executing {cmd=}")
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    initial_dbg_path = asm_path / cfg.full_asm_subdir / "00_CoverageBasedCorrection" / "initial_dbg.gfa"
    assert initial_dbg_path.exists(), f"Initial DBG file {initial_dbg_path} does not exist"
    mapped_ids = map_lja_ids_to_jumbo_ids(initial_dbg_path, graph_dbg_path)

    for stage in cfg.eval_stages:
        logger.info(f"Evaluation {stage=}...")
        alignments_path = asm_path / stage / "alignments.txt"
        assert alignments_path.exists(), f"Alignments file {alignments_path} does not exist"

        evaluation = evaluate_la_jolla(mult_info_path, alignments_path, mapped_ids)

        output_path = output_path_eval / stage
        if not output_path.exists():
            output_path.mkdir(parents=True)

        evaluation_file = output_path / "evaluation.json"
        logger.info(f"Writing evaluation data to {evaluation_file}.")
        with open(evaluation_file, "w") as f:
            json.dump(evaluation, f, indent=4)


def config_check(cfg: DictConfig):
    assert cfg.threads is not None, "Number of threads is not defined!"
    assert Path(cfg.eval_path).exists(), f"{cfg.eval_path} does not exist!"
    assert Path(cfg.lja_bin_path).exists(), f"{cfg.lja_bin_path} does not exist!"


@hydra.main(version_base=None, config_path=str(ph.get_config_root()), config_name="eval.yaml")
def main(cfg: DictConfig):
    config_check(cfg)

    for subdir in os.listdir(cfg.eval_path):
        try:
            # This works for positive integers, but that is our only use case so far
            if not subdir.isdigit():
                continue
            eval_lja(cfg, subdir=subdir)
        except Exception as ex:
            logger.warning(f"Exception {ex}")
            logger.warning(f"Failed evaluating {subdir} sample, continue with next sample...")

    eval_path = Path(cfg.eval_path)
    df = get_eval_table(eval_path)
    df.to_csv(eval_path / "eval_summary.csv", index=None)


if __name__ == "__main__":
    main()
