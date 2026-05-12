"""Tests for the kNN module — format-level only.

scanpy's neighbor computation itself is covered upstream; here we only
exercise the on-disk Neighbors format (shared cell_ids + two CSR groups).
"""

import h5py
import numpy as np
import scipy.sparse as sp

from schemas import Neighbors


def _toy_neighbors():
    ids = ["a", "b", "c", "d"]
    d = sp.csr_matrix(np.array([[0, 1, 2, 0], [1, 0, 0, 3], [2, 0, 0, 4], [0, 3, 4, 0]], dtype=float))
    c = sp.csr_matrix(np.array([[0, 0.5, 0.2, 0], [0.5, 0, 0, 0.1], [0.2, 0, 0, 0.9], [0, 0.1, 0.9, 0]]))
    return Neighbors(cell_ids=ids, distances=d, connectivities=c)


def test_neighbors_roundtrip(tmp_path):
    nbrs = _toy_neighbors()
    p = tmp_path / "neighbors.h5"
    nbrs.write(p)
    loaded = Neighbors.read(p)

    assert loaded.cell_ids == nbrs.cell_ids
    np.testing.assert_array_equal(loaded.distances.toarray(), nbrs.distances.toarray())
    np.testing.assert_array_equal(
        loaded.connectivities.toarray(), nbrs.connectivities.toarray()
    )


def test_neighbors_h5_layout(tmp_path):
    """cell_ids must be stored once at the top level, with both CSR groups beside it."""
    nbrs = _toy_neighbors()
    p = tmp_path / "neighbors.h5"
    nbrs.write(p)

    with h5py.File(p, "r") as h5:
        assert set(h5.keys()) == {"cell_ids", "distances", "connectivities"}
        for group in ("distances", "connectivities"):
            assert set(h5[group].keys()) == {"data", "indices", "indptr"}
