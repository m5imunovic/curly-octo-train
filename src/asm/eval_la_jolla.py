"""Evaluates the corrected reads producesd by La Jolla Assembler. More specifically, it evaluates
the corrected reads produced by the second step of the pipeline (i.e. the topology-based correction).
"""
import os
import json
from pathlib import Path
from typing import Dict

import hydra
from omegaconf import DictConfig
from typeguard import typechecked

import utils.path_helpers as ph
from asm.lja_rc_map import get_rc_map_mp_pool_batch
from asm.mult_info_parser import parse_mult_info, partition_mult_info_edges



def parse_alignments_entry(read_id: str, edge_ids: str):
    """First line is read ID, second line is list of edge IDs. 
    Returns a tuple of (read_id, edge_ids, is_rc) where is_rc is a list of 0s and 1s indicating
    whether the corresponding edge string is reverse complemented."""
    read_id = read_id.strip()[1:]
    edge_ids = edge_ids.strip().split()
    is_rc = [edge_id.startswith('-') for edge_id in edge_ids]
    edge_ids = [edge_id[1:] if is_rc else edge_id for edge_id, is_rc in zip(edge_ids, is_rc)]

    return read_id, edge_ids, is_rc
           

@typechecked
def parse_alignments(alignments_path: Path) -> set:
    with open(alignments_path) as f:
        correct_reads = set()
        for read_line in f:
            _, edge_ids, _ = parse_alignments_entry(read_line, next(f))
            correct_reads.update(edge_ids)

    return correct_reads


@typechecked
def get_correct_edges(alignments_path: Path, rc_map: Dict) -> set:
    correct_edges = parse_alignments(alignments_path)
    # if read_id does not exist we simply re-add it as we are doing union with input set anyway
    correct_edges_rc = set(rc_map.get(read_id, read_id) for read_id in correct_edges)

    return correct_edges | correct_edges_rc


@typechecked
def get_confusion_matrix(mult_info_path: Path, alignments_path: Path, rc_map: Dict) -> tuple:
    """
    Requires the following files:
    - alignments.txt
    - mult.info
    - initial_dbg.gfa
    """

    # get ground truth
    mult_info = parse_mult_info(mult_info_path)
    correct_edges_gt, incorrect_edges_gt = partition_mult_info_edges(mult_info)
    print("Ground truth:")
    print(f"Correct edges: {len(correct_edges_gt)}")
    print(f"Incorrect edges: {len(incorrect_edges_gt)}")

    # get mowerDBG edge assignments
    correct_edges = get_correct_edges(alignments_path, rc_map=rc_map)
    incorrect_edges = (correct_edges_gt | incorrect_edges_gt) - correct_edges
    print("MowerDBG classified:")
    print(f"\t{len(correct_edges)} as correct edges.")
    print(f"\t{len(incorrect_edges)} as incorrect egdes.")

    # true positives
    tp = len(correct_edges_gt & correct_edges)
    tn = len(incorrect_edges_gt & incorrect_edges)
    fp = len(incorrect_edges_gt & correct_edges)
    fn = len(correct_edges_gt & incorrect_edges)

    return tp, tn, fp, fn


@typechecked
def calculate_metrics(tp: int, tn:int, fp:int , fn:int) -> Dict[str, float]:
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


def evaluate_la_jolla(mult_info_path: Path, alignments_path: Path, gfa_path: Path, k: int, threads: int) -> Dict:

    rc_map = get_rc_map_mp_pool_batch(gfa_path=gfa_path, k=k, threads=threads)
    tp, tn, fp, fn = get_confusion_matrix(mult_info_path, alignments_path=alignments_path, rc_map=rc_map)
    metrics = calculate_metrics(tp, tn, fp, fn)

    evaluations = {
        "Mult info path": str(mult_info_path),
        "Alignments path": str(alignments_path),
        "Confusion matrix": {
            "tp fn": [tp, fn],
            "fp tn": [fp, tn],
        },
        "metrics": metrics
    }

    return evaluations



@hydra.main(version_base="1.2", config_path=ph.get_config_root(), config_name="eval.yaml")
def main(cfg: DictConfig):
    mult_info_path = Path(cfg.asm_path) / "mult.info"
    assert mult_info_path.exists(), f"Mult info file {mult_info_path} does not exist"
    alignments_path = Path(cfg.asm_path) / cfg.eval_stage / "alignments.txt"
    assert alignments_path.exists(), f"Alignments file {alignments_path} does not exist"
    initial_dbg_path = Path(cfg.asm_path) / cfg.full_asm_subdir / "00_CoverageBasedCorrection" / "initial_dbg.gfa"
    assert initial_dbg_path.exists(), f"Initial DBG file {initial_dbg_path} does not exist"
    k = cfg.k
    threads = cfg.threads or os.cpu_count() - 1
    output_path = Path(cfg.output_path)
    evaluation = evaluate_la_jolla(mult_info_path, alignments_path, gfa_path=initial_dbg_path, k=k, threads=threads)
    if not output_path.exists():
        output_path.mkdir(parents=True)

    with open(output_path / "evaluation.json", "w") as f:
        json.dump(evaluation, f, indent=4)


if __name__ == "__main__":
    main()
