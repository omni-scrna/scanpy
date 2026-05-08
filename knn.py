#!/usr/bin/env python3
"""kNN graph module (scanpy-backed) for omnibenchmark.

Outputs two HDF5 files per run:
  {name}_distances.h5      — KNN distance matrix (sparse CSR)
  {name}_connectivities.h5 — UMAP connectivities matrix (sparse CSR)

Each file layout:
  /cell_ids  string array (n_cells,)
  /data      float64 array (nnz,)
  /indices   int32 array (nnz,)
  /indptr    int32 array (n_cells+1,)
"""

import sys
from pathlib import Path

import scipy.sparse as sp
import anndata as ad
import h5py
import numpy as np
import polars as pl
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_knn_parser  # noqa: E402


def write_sparse(h5, m: sp.spmatrix, cell_ids: list[str]) -> None:
    m = m.tocsr()
    h5.create_dataset("cell_ids", data=np.array(cell_ids, dtype=bytes))
    h5.create_dataset("data", data=m.data)
    h5.create_dataset("indices", data=m.indices)
    h5.create_dataset("indptr", data=m.indptr)


def main():
    args = build_knn_parser().parse_args()

    df = pl.read_csv(args.pcas, separator="\t")
    cell_ids = df["cell_id"].to_list()
    embedding = df.drop("cell_id").to_numpy().astype(np.float64)

    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
    adata.obs_names = cell_ids
    adata.obsm["X_pca"] = embedding

    sc.pp.neighbors(
        adata,
        n_neighbors=args.n_neighbors,
        method=args.flavor,
        use_rep="X_pca",
        random_state=args.random_seed,
    )

    for key in ("distances", "connectivities"):
        out = Path(args.output_dir) / f"{args.name}_{key}.h5"
        with h5py.File(out, "w") as h5:
            write_sparse(h5, adata.obsp[key], cell_ids)


if __name__ == "__main__":
    main()
