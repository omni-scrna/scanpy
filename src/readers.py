"""Composite readers — pure deserialisation lives on the schemas themselves."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from schemas import Neighbors


def read_neighbors_as_anndata(
    path: str | Path, format: str = "h5"
) -> tuple[ad.AnnData, list[str]]:
    """Load a Neighbors bundle and wrap it in an AnnData ready for scanpy graph algos."""
    nbrs = Neighbors.read(path, format=format)
    adata = ad.AnnData(obs=pd.DataFrame(index=nbrs.cell_ids))
    adata.obsp["distances"] = nbrs.distances
    adata.obsp["connectivities"] = nbrs.connectivities
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    return adata, nbrs.cell_ids
