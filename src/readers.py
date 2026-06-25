"""Readers for artifacts produced by this module, conforming to the omni-scrna spec."""

from pathlib import Path

import h5py
import scipy.sparse as sp


def _read_csr(grp, shape):
    return sp.csr_matrix(
        (grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape
    )


def read_neighbors(path):
    """Load {name}_neighbors.h5 written by knn.py: the kNN distance graph (flat CSR
    at the file root) and the connectivities (nested /connectivities group).

    Returns (distances, connectivities, cell_ids): two square CSR matrices sharing
    the returned cell ids.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        cell_ids = [c.decode() if isinstance(c, bytes) else c for c in h5["cell_ids"][:]]
        n = len(cell_ids)
        distances = _read_csr(h5, (n, n))
        connectivities = _read_csr(h5["connectivities"], (n, n))
    return distances, connectivities, cell_ids
