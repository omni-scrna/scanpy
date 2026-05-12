"""Composite readers — pure deserialisation lives on the schemas themselves."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
import pandas as pd
import scipy.sparse as sp

from schemas import Neighbors


def read_tenx_h5(path: str | Path) -> ad.AnnData:
    """Load a TENx-format HDF5 matrix (genes x cells, CSC) as a cells x genes AnnData.

    Expected layout under ``/matrix``: ``data``, ``indices``, ``indptr``,
    ``shape``, ``genes``, ``barcodes``.
    """
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])
        gene_ids = g["genes"][:].astype(str)
        cell_ids = g["barcodes"][:].astype(str)

    X = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    adata = ad.AnnData(X=X)
    adata.obs_names = cell_ids
    adata.var_names = gene_ids
    return adata


def read_neighbors_as_anndata(
    path: str | Path,
) -> tuple[ad.AnnData, list[str]]:
    """Load a Neighbors bundle and wrap it in an AnnData ready for scanpy graph algos."""
    nbrs = Neighbors.read(path)
    adata = ad.AnnData(obs=pd.DataFrame(index=pd.Index(nbrs.cell_ids)))
    adata.obsp["distances"] = nbrs.distances
    adata.obsp["connectivities"] = nbrs.connectivities
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    return adata, nbrs.cell_ids
