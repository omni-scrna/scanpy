#!/usr/bin/env python3
"""kNN graph module (scanpy-backed) for omnibenchmark.

Outputs a single HDF5 file per run:
  {output_dir}/{name}_neighbors.h5 — KNN distances + UMAP connectivities sharing one cell_ids list.

File layout:
  /cell_ids                  string array (n_cells,)
  /distances/{data,indices,indptr}       sparse CSR (n_cells, n_cells)
  /connectivities/{data,indices,indptr}  sparse CSR (n_cells, n_cells)
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import polars as pl
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_knn_parser  # noqa: E402
from writers import Neighbors, write_neighbors  # noqa: E402


def compute_neighbors(embedding, cell_ids, n_neighbors, flavor, random_seed):
    """Run scanpy's neighbors on an embedding and return a Neighbors object."""
    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
    adata.obs_names = cell_ids
    adata.obsm["X_pca"] = embedding

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        method=flavor,
        use_rep="X_pca",
        random_state=random_seed,
    )
    return Neighbors(
        cell_ids=cell_ids,
        distances=adata.obsp["distances"],
        connectivities=adata.obsp["connectivities"],
    )


def main():
    args = build_knn_parser().parse_args()

    df = pl.read_csv(args.pcas, separator="\t")
    cell_ids = df["cell_id"].to_list()
    embedding = df.drop("cell_id").to_numpy().astype(np.float64)

    nbrs = compute_neighbors(
        embedding,
        cell_ids,
        n_neighbors=args.n_neighbors,
        flavor=args.flavor,
        random_seed=args.random_seed,
    )

    out = Path(args.output_dir) / f"{args.name}_neighbors.h5"
    write_neighbors(nbrs, out)


if __name__ == "__main__":
    main()
