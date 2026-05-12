"""Reusable readers and writers for objects produced by this module, conforming to the omni-scrna benchmark spec."""

from dataclasses import dataclass, field

import h5py
import numpy as np
import polars as pl
import scipy.sparse as sp


@dataclass
class Embedding:
    matrix: np.ndarray  # shape (n_cells, n_dims)
    row_ids: list  # cell barcodes, length n_cells
    col_names: list = field(default_factory=list)  # dim labels; auto-generated if empty


@dataclass
class SparseGraph:
    matrix: sp.csr_matrix  # square (n_cells, n_cells)
    cell_ids: list[str]


@dataclass
class Neighbors:
    cell_ids: list[str]
    distances: sp.csr_matrix  # square (n_cells, n_cells)
    connectivities: sp.csr_matrix  # square (n_cells, n_cells)


@dataclass
class Clustering:
    cell_ids: list[str]
    labels: list[str]


# ── Embedding ──────────────────────────────────────────────────────────────────


def _col_names(embedding):
    if embedding.col_names:
        return embedding.col_names
    return [f"dim_{i + 1}" for i in range(embedding.matrix.shape[1])]


def _write_embedding_tsv(path, embedding):
    cols = _col_names(embedding)
    df = pl.from_numpy(embedding.matrix, schema=cols).insert_column(
        0, pl.Series("", embedding.row_ids)
    )
    with open(path, "w") as f:
        f.write("cell_id\t" + "\t".join(cols) + "\n")
        df.write_csv(f, separator="\t", include_header=False)


def write_embeddings(obj, path, format="tsv"):
    """Write an Embedding.

    tsv layout:
      header: cell_id<TAB>dim_1<TAB>...<TAB>dim_n   (col_names override dim_*)
      body:   one row per cell, tab-separated.
    """
    if format == "tsv":
        _write_embedding_tsv(path, obj)
    else:
        raise ValueError(f"unsupported format: {format!r}")


# ── SparseGraph ────────────────────────────────────────────────────────────────


def _write_graph_h5(path, graph):
    m = graph.matrix.tocsr()
    with h5py.File(path, "w") as h5:
        h5.create_dataset("cell_ids", data=np.array(graph.cell_ids, dtype=bytes))
        h5.create_dataset("data", data=m.data)
        h5.create_dataset("indices", data=m.indices)
        h5.create_dataset("indptr", data=m.indptr)


def _read_graph_h5(path):
    with h5py.File(path, "r") as h5:
        cell_ids = [
            s.decode() if isinstance(s, bytes) else s for s in h5["cell_ids"][:]
        ]
        mat = sp.csr_matrix(
            (h5["data"][:], h5["indices"][:], h5["indptr"][:]),
            shape=(len(cell_ids), len(cell_ids)),
        )
    return SparseGraph(matrix=mat, cell_ids=cell_ids)


def write_graph(obj, path, format="h5"):
    """Write a SparseGraph.

    h5 layout:
      /cell_ids   string array (n,)
      /data       float array (nnz,)
      /indices    int array   (nnz,)
      /indptr     int array   (n+1,)
      CSR matrix is square (n, n).
    """
    if format == "h5":
        _write_graph_h5(path, obj)
    else:
        raise ValueError(f"unsupported format: {format!r}")


def read_graph(path, format="h5"):
    """Inverse of :func:`write_graph`. See it for the on-disk layout."""
    if format == "h5":
        return _read_graph_h5(path)
    else:
        raise ValueError(f"unsupported format: {format!r}")


# ── Neighbors ──────────────────────────────────────────────────────────────────


def _write_csr_group(group, m):
    m = m.tocsr()
    group.create_dataset("data", data=m.data)
    group.create_dataset("indices", data=m.indices)
    group.create_dataset("indptr", data=m.indptr)


def _read_csr_group(group, n):
    return sp.csr_matrix(
        (group["data"][:], group["indices"][:], group["indptr"][:]),
        shape=(n, n),
    )


def _write_neighbors_h5(path, neighbors):
    with h5py.File(path, "w") as h5:
        h5.create_dataset("cell_ids", data=np.array(neighbors.cell_ids, dtype=bytes))
        _write_csr_group(h5.create_group("distances"), neighbors.distances)
        _write_csr_group(h5.create_group("connectivities"), neighbors.connectivities)


def _read_neighbors_h5(path):
    with h5py.File(path, "r") as h5:
        cell_ids = [
            s.decode() if isinstance(s, bytes) else s for s in h5["cell_ids"][:]
        ]
        n = len(cell_ids)
        return Neighbors(
            cell_ids=cell_ids,
            distances=_read_csr_group(h5["distances"], n),
            connectivities=_read_csr_group(h5["connectivities"], n),
        )


def write_neighbors(obj, path, format="h5"):
    """Write a Neighbors bundle.

    h5 layout:
      /cell_ids                              string array (n,)
      /distances/{data,indices,indptr}       sparse CSR (n, n)
      /connectivities/{data,indices,indptr}  sparse CSR (n, n)
    Both matrices share the top-level /cell_ids.
    """
    if format == "h5":
        _write_neighbors_h5(path, obj)
    else:
        raise ValueError(f"unsupported format: {format!r}")


def read_neighbors(path, format="h5"):
    """Inverse of :func:`write_neighbors`. See it for the on-disk layout."""
    if format == "h5":
        return _read_neighbors_h5(path)
    else:
        raise ValueError(f"unsupported format: {format!r}")


# ── Clustering ─────────────────────────────────────────────────────────────────


def _write_clustering_tsv(path, clustering):
    pl.DataFrame(
        {"cell_id": clustering.cell_ids, "cluster": list(clustering.labels)}
    ).write_csv(path, separator="\t")


def write_clustering(obj, path, format="tsv"):
    """Write a Clustering.

    tsv layout:
      header: cell_id<TAB>cluster
      body:   one row per cell.
    """
    if format == "tsv":
        _write_clustering_tsv(path, obj)
    else:
        raise ValueError(f"unsupported format: {format!r}")
