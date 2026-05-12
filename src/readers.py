"""Reusable readers that hydrate scanpy/anndata objects from on-disk artefacts."""

import anndata as ad
import pandas as pd

from writers import read_neighbors


def read_neighbors_as_anndata(path, format="h5"):
    """Load a Neighbors bundle and wrap it in an AnnData ready for scanpy graph algos."""
    nbrs = read_neighbors(path, format=format)
    adata = ad.AnnData(obs=pd.DataFrame(index=nbrs.cell_ids))
    adata.obsp["distances"] = nbrs.distances
    adata.obsp["connectivities"] = nbrs.connectivities
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    return adata, nbrs.cell_ids
