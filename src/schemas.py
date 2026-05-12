"""Data schemas exchanged between omni-scrna scanpy stages.

Each dataclass owns its on-disk layout via ``read`` / ``write`` methods.
Layouts are documented on the methods themselves.

TODO: if/when a second format is added per schema, reintroduce a
``format=`` parameter (or a small dispatch mixin) — see git history.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
import polars as pl
import scipy.sparse as sp

PathLike = str | Path


# ── helpers ───────────────────────────────────────────────────────────────────


def _write_csr_group(group: h5py.Group, m: sp.csr_matrix) -> None:
    m = m.tocsr()
    group.create_dataset("data", data=m.data)
    group.create_dataset("indices", data=m.indices)
    group.create_dataset("indptr", data=m.indptr)


def _read_csr_group(group: h5py.Group, n: int) -> sp.csr_matrix:
    return sp.csr_matrix(
        (group["data"][:], group["indices"][:], group["indptr"][:]),
        shape=(n, n),
    )


def _decode_cell_ids(dataset: h5py.Dataset) -> list[str]:
    return [s.decode() if isinstance(s, bytes) else s for s in dataset[:]]


# ── Embedding ──────────────────────────────────────────────────────────────────


@dataclass
class Embedding:
    matrix: np.ndarray  # shape (n_cells, n_dims)
    row_ids: list[str]  # cell barcodes, length n_cells
    col_names: list[str] = field(default_factory=list)  # dim labels; auto if empty

    def write(self, path: PathLike) -> None:
        """tsv layout: header ``cell_id<TAB>dim_1<TAB>...``; one row per cell."""
        cols = self.col_names or [f"dim_{i + 1}" for i in range(self.matrix.shape[1])]
        df = pl.from_numpy(self.matrix, schema=cols).insert_column(
            0, pl.Series("", self.row_ids)
        )
        with open(path, "w") as f:
            f.write("cell_id\t" + "\t".join(cols) + "\n")
            df.write_csv(f, separator="\t", include_header=False)

    @classmethod
    def read(cls, path: PathLike) -> Embedding:
        df = pl.read_csv(path, separator="\t")
        return cls(
            matrix=df.drop("cell_id").to_numpy().astype(np.float64),
            row_ids=df["cell_id"].to_list(),
            col_names=[c for c in df.columns if c != "cell_id"],
        )


# ── SparseGraph ────────────────────────────────────────────────────────────────


@dataclass
class SparseGraph:
    matrix: sp.csr_matrix  # square (n_cells, n_cells)
    cell_ids: list[str]

    def write(self, path: PathLike) -> None:
        """h5 layout: /cell_ids + flat CSR (/data, /indices, /indptr)."""
        with h5py.File(path, "w") as h5:
            h5.create_dataset("cell_ids", data=np.array(self.cell_ids, dtype=bytes))
            _write_csr_group(h5, self.matrix)

    @classmethod
    def read(cls, path: PathLike) -> SparseGraph:
        with h5py.File(path, "r") as h5:
            cell_ids = _decode_cell_ids(h5["cell_ids"])
            mat = _read_csr_group(h5, len(cell_ids))
        return cls(matrix=mat, cell_ids=cell_ids)


# ── Neighbors ──────────────────────────────────────────────────────────────────


@dataclass
class Neighbors:
    cell_ids: list[str]
    distances: sp.csr_matrix  # square (n_cells, n_cells)
    connectivities: sp.csr_matrix  # square (n_cells, n_cells)

    def write(self, path: PathLike) -> None:
        """h5 layout: /cell_ids + /distances/{data,indices,indptr} + /connectivities/{...}."""
        with h5py.File(path, "w") as h5:
            h5.create_dataset("cell_ids", data=np.array(self.cell_ids, dtype=bytes))
            _write_csr_group(h5.create_group("distances"), self.distances)
            _write_csr_group(h5.create_group("connectivities"), self.connectivities)

    @classmethod
    def read(cls, path: PathLike) -> Neighbors:
        with h5py.File(path, "r") as h5:
            cell_ids = _decode_cell_ids(h5["cell_ids"])
            n = len(cell_ids)
            return cls(
                cell_ids=cell_ids,
                distances=_read_csr_group(h5["distances"], n),
                connectivities=_read_csr_group(h5["connectivities"], n),
            )


# ── Clustering ─────────────────────────────────────────────────────────────────


@dataclass
class Clustering:
    cell_ids: list[str]
    labels: list[str]

    def write(self, path: PathLike) -> None:
        """tsv layout: header ``cell_id<TAB>cluster``; one row per cell."""
        pl.DataFrame(
            {"cell_id": self.cell_ids, "cluster": list(self.labels)}
        ).write_csv(path, separator="\t")

    @classmethod
    def read(cls, path: PathLike) -> Clustering:
        df = pl.read_csv(path, separator="\t")
        return cls(
            cell_ids=df["cell_id"].to_list(),
            labels=df["cluster"].cast(pl.Utf8).to_list(),
        )
