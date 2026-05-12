"""Tests for the clustering module."""

import numpy as np
import polars as pl
import pytest
import scipy.sparse as sp

from cluster import cluster_leiden
from readers import read_graph_as_anndata
from writers import SparseGraph, write_clustering, write_graph, Clustering


def _two_blobs_graph(n_per_block=8):
    """Build a connectivities matrix with two well-separated cliques."""
    n = 2 * n_per_block
    a = np.ones((n_per_block, n_per_block))
    block = np.block([[a, np.zeros_like(a)], [np.zeros_like(a), a]])
    np.fill_diagonal(block, 0)
    ids = [f"c{i}" for i in range(n)]
    return SparseGraph(matrix=sp.csr_matrix(block), cell_ids=ids)


@pytest.fixture
def graph_h5(tmp_path):
    g = _two_blobs_graph()
    p = tmp_path / "graph.h5"
    write_graph(g, p)
    return p, g


def test_read_graph_as_anndata_shape(graph_h5):
    path, g = graph_h5
    adata, cell_ids = read_graph_as_anndata(path)
    assert cell_ids == g.cell_ids
    assert adata.n_obs == len(g.cell_ids)
    assert "connectivities" in adata.obsp
    assert adata.uns["neighbors"]["connectivities_key"] == "connectivities"


def test_cluster_leiden_separates_blocks(graph_h5):
    path, g = graph_h5
    adata, _ = read_graph_as_anndata(path)
    labels = cluster_leiden(adata, resolution=1.0, random_seed=0)

    assert len(labels) == len(g.cell_ids)
    n = len(g.cell_ids) // 2
    block_a = set(labels[:n])
    block_b = set(labels[n:])
    assert len(block_a) == 1 and len(block_b) == 1
    assert block_a != block_b


def test_cluster_leiden_deterministic(graph_h5):
    path, _ = graph_h5
    adata1, _ = read_graph_as_anndata(path)
    adata2, _ = read_graph_as_anndata(path)
    l1 = list(cluster_leiden(adata1, random_seed=42))
    l2 = list(cluster_leiden(adata2, random_seed=42))
    assert l1 == l2


def test_write_clustering_tsv(tmp_path):
    out = tmp_path / "clusters.tsv"
    write_clustering(Clustering(cell_ids=["a", "b", "c"], labels=["0", "0", "1"]), out)
    df = pl.read_csv(out, separator="\t")
    assert df.columns == ["cell_id", "cluster"]
    assert df["cell_id"].to_list() == ["a", "b", "c"]
    assert df["cluster"].cast(pl.Utf8).to_list() == ["0", "0", "1"]
