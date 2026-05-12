"""Tests for the clustering module."""

import numpy as np
import polars as pl
import pytest
import scipy.sparse as sp

from cluster import cluster_leiden
from readers import read_neighbors_as_anndata
from schemas import Clustering, Neighbors


def _two_blobs_csr(n_per_block=8):
    n = 2 * n_per_block
    a = np.ones((n_per_block, n_per_block))
    block = np.block([[a, np.zeros_like(a)], [np.zeros_like(a), a]])
    np.fill_diagonal(block, 0)
    return sp.csr_matrix(block), [f"c{i}" for i in range(n)]


@pytest.fixture
def neighbors_h5(tmp_path):
    mat, ids = _two_blobs_csr()
    nbrs = Neighbors(cell_ids=ids, distances=mat.copy(), connectivities=mat.copy())
    p = tmp_path / "neighbors.h5"
    nbrs.write(p)
    return p, nbrs


def test_read_neighbors_as_anndata_shape(neighbors_h5):
    path, nbrs = neighbors_h5
    adata, cell_ids = read_neighbors_as_anndata(path)
    assert cell_ids == nbrs.cell_ids
    assert adata.n_obs == len(nbrs.cell_ids)
    assert "distances" in adata.obsp
    assert "connectivities" in adata.obsp
    assert adata.uns["neighbors"]["connectivities_key"] == "connectivities"


def test_cluster_leiden_separates_blocks(neighbors_h5):
    path, nbrs = neighbors_h5
    adata, _ = read_neighbors_as_anndata(path)
    labels = cluster_leiden(adata, resolution=1.0, random_seed=0)

    assert len(labels) == len(nbrs.cell_ids)
    n = len(nbrs.cell_ids) // 2
    block_a = set(labels[:n])
    block_b = set(labels[n:])
    assert len(block_a) == 1 and len(block_b) == 1
    assert block_a != block_b


def test_cluster_leiden_deterministic(neighbors_h5):
    path, _ = neighbors_h5
    adata1, _ = read_neighbors_as_anndata(path)
    adata2, _ = read_neighbors_as_anndata(path)
    l1 = list(cluster_leiden(adata1, random_seed=42))
    l2 = list(cluster_leiden(adata2, random_seed=42))
    assert l1 == l2


def test_write_clustering_tsv(tmp_path):
    out = tmp_path / "clusters.tsv"
    Clustering(cell_ids=["a", "b", "c"], labels=["0", "0", "1"]).write(out)
    df = pl.read_csv(out, separator="\t")
    assert df.columns == ["cell_id", "cluster"]
    assert df["cell_id"].to_list() == ["a", "b", "c"]
    assert df["cluster"].cast(pl.Utf8).to_list() == ["0", "0", "1"]
