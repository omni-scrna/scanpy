"""Clustering reads the stored distances + connectivities graph and Leiden
recovers the structure that is actually in the data.

The graph is built from the be1 fixture's real PCA (conftest.be1_pcas_tsv), not
from gaussian blobs pushed 25 sigma apart: on synthetic blobs "Leiden separates
them" is true for any graph that is not outright corrupt, so the test could not
fail for any interesting reason. On the eight real cell lines it can.
"""

import anndata as ad
import numpy as np
import polars as pl
import pytest
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from cluster import build_adata, cluster_leiden
from knn import write_neighbors_graph

K = 15
LEIDEN = ("igraph", "RBConfiguration", 1.0, 0)


@pytest.fixture(scope="module")
def neighbors_h5(be1_pcas_tsv, tmp_path_factory):
    df = pl.read_csv(be1_pcas_tsv, separator="\t")
    ids = df["cell_id"].to_list()

    a = ad.AnnData(X=np.zeros((len(ids), 1)))
    a.obs_names = ids
    a.obsm["X_pca"] = df.drop("cell_id").to_numpy()
    sc.pp.neighbors(a, n_neighbors=K, method="umap", use_rep="X_pca", random_state=0)

    tmp = tmp_path_factory.mktemp("clust")
    write_neighbors_graph(a, tmp, "t")
    return tmp / "t_neighbors.h5", ids


def test_cluster_leiden_recovers_the_cell_lines(neighbors_h5, be1_truth):
    path, ids = neighbors_h5
    adata, cell_ids = build_adata(path)
    labels = cluster_leiden(adata, *LEIDEN)
    assert len(labels) == len(ids)

    truth = [be1_truth[c] for c in cell_ids]
    # Leiden at resolution 1.0 splits some cell lines in two, which is a
    # resolution choice, not a failure. What must never happen is two cell lines
    # collapsing into one cluster, so pin the dominant cluster per label.
    modal = {}
    for lab in set(truth):
        counts = np.bincount([int(c) for t, c in zip(truth, labels) if t == lab])
        modal[lab] = counts.argmax()
    assert len(set(modal.values())) == len(modal), f"cell lines share a cluster: {modal}"
    # ...and a floor on the partition as a whole, so shredding it into fragments
    # (which would also give every label its own modal cluster) still fails.
    assert adjusted_rand_score(truth, labels) > 0.6


def test_cluster_leiden_deterministic(neighbors_h5):
    path, _ = neighbors_h5
    a1, _ = build_adata(path)
    a2, _ = build_adata(path)
    assert cluster_leiden(a1, *LEIDEN) == cluster_leiden(a2, *LEIDEN)
