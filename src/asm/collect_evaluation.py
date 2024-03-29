"""Collects (LJA) evaluation results from eval folder."""
import json
import sys
from pathlib import Path

import pandas as pd


def get_eval_table(eval_path: Path):
    """Finds all `evaluation.json` files under `eval_path` directory and extracts the data, adding one row per file
    into the data table.

    Expected columns are:
    'sample', 'phase', 'Precision', 'Recall', 'Accuracy', 'F1', 'tp', 'fn', 'fp', 'tn', 'edge_cnt'
    """
    eval_files = eval_path.glob("**/evaluation.json")

    data = []
    columns = ["sample", "phase"]
    columns_ready = False
    for eval_file in sorted(eval_files):
        values = []
        values.append(eval_file.parent.parent.parent.name)  # sample name
        values.append(eval_file.parent.name)  # phase name

        with open(eval_file) as f:
            evaluation_data = json.load(f)
            confusion_matrix = evaluation_data["Confusion matrix"]

            metrics = evaluation_data["metrics"]
            for col, val in metrics.items():
                if not columns_ready:
                    columns.append(col)
                values.append(f"{val:.04f}")

            edge_cnt = 0
            for col, val in confusion_matrix.items():
                key1, key2 = col.split(" ")
                if not columns_ready:
                    columns.append(key1)
                    columns.append(key2)
                values.append(val[0])
                values.append(val[1])
                edge_cnt += sum(val)
            if not columns_ready:
                columns.append("edge_cnt")
            values.append(edge_cnt)

        data.append(values)
        columns_ready = True

    df = pd.DataFrame(data=data, columns=columns, index=None)
    return df


if __name__ == "__main__":
    # TODO: call this at the end of evaluation loop instead
    eval_path = Path(sys.argv[1])
    assert eval_path.exists(), f"{eval_path} does not exists"
    df = get_eval_table(eval_path)
    df.to_csv(eval_path / "eval_summary.csv", index=None)
