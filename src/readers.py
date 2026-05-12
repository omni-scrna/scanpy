"""Reusable readers that hydrate scanpy/anndata objects from on-disk artefacts."""

import anndata as ad
import pandas as pd

from writers import read_graph


def read_graph_as_anndata(path, format="h5"):
    """Load a SparseGraph and wrap it in an AnnData ready for scanpy graph algos."""
    graph = read_graph(path, format=format)
    adata = ad.AnnData(obs=pd.DataFrame(index=graph.cell_ids))
    adata.obsp["connectivities"] = graph.matrix
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}
    return adata, graph.cell_ids
