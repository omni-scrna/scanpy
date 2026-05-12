"""Tests for the PCA module — format-level + minimal wiring.

scanpy's PCA itself is upstream; we exercise:
  - Embedding TSV round-trip (the on-disk contract)
  - run_pca returns an Embedding with the right shape / cell_ids
"""

import argparse

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp

from pca import run_pca
from schemas import Embedding


def _toy_adata(n_cells: int = 30, n_genes: int = 20, seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.random((n_cells, n_genes)))
    a = ad.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(n_cells)]
    a.var_names = [f"g{i}" for i in range(n_genes)]
    return a


def test_embedding_tsv_roundtrip(tmp_path):
    emb = Embedding(
        matrix=np.array([[1.0, 2.0], [3.0, 4.0]]),
        row_ids=["a", "b"],
        col_names=["PC1", "PC2"],
    )
    p = tmp_path / "emb.tsv"
    emb.write(p)
    loaded = Embedding.read(p)
    np.testing.assert_array_equal(loaded.matrix, emb.matrix)
    assert loaded.row_ids == emb.row_ids
    assert loaded.col_names == emb.col_names


def test_embedding_tsv_default_dim_names(tmp_path):
    """Empty col_names => header uses auto dim_1..dim_n."""
    emb = Embedding(matrix=np.zeros((2, 3)), row_ids=["a", "b"], col_names=[])
    p = tmp_path / "emb.tsv"
    emb.write(p)
    df = pl.read_csv(p, separator="\t")
    assert df.columns == ["cell_id", "dim_1", "dim_2", "dim_3"]


def test_run_pca_returns_embedding_with_expected_shape():
    adata = _toy_adata(n_cells=30, n_genes=20)
    args = argparse.Namespace(n_components=5, solver="arpack", random_seed=0)
    emb = run_pca(adata, args)
    assert isinstance(emb, Embedding)
    assert emb.matrix.shape == (30, 5)
    assert emb.row_ids == list(adata.obs_names)
    assert emb.col_names == ["PC1", "PC2", "PC3", "PC4", "PC5"]
