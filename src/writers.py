"""Reusable writers for objects produced by this module, conforming to the omni-scrna benchmark spec."""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

@dataclass
class Embedding:
    matrix: np.ndarray        # shape (n_cells, n_dims)
    row_ids: list             # cell barcodes, length n_cells
    col_names: list = field(default_factory=list)  # dim labels; auto-generated if empty


def _col_names(embedding):
    if embedding.col_names:
        return embedding.col_names
    return [f"dim_{i + 1}" for i in range(embedding.matrix.shape[1])]


def _write_tsv(path, embedding):
    # index=True + no index name → header row has N cols, data rows N+1,
    # which read.table(f, header=TRUE) auto-promotes to row.names.
    df = pd.DataFrame(embedding.matrix, index=embedding.row_ids, columns=_col_names(embedding))
    df.to_csv(path, sep="\t")


def write_embeddings(obj, path, format="tsv"):
    if format == "tsv":
        _write_tsv(path, obj)
    else:
        raise ValueError(f"unsupported format: {format!r}")
