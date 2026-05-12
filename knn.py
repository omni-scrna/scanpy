#!/usr/bin/env python3
"""kNN graph module (scanpy-backed) for omnibenchmark.

Input:
  --pcas.tsv — PCA embedding TSV (produced by ``pca.py``)
    header: cell_id<TAB>dim_1<TAB>...<TAB>dim_n
    body:   one row per cell; numeric columns parsed as float64.

Outputs a single HDF5 file per run:
  {output_dir}/{name}_neighbors.h5 — KNN distances + UMAP connectivities sharing one cell_ids list.

File layout:
  /cell_ids                              string array (n_cells,)
  /distances/{data,indices,indptr}       sparse CSR (n_cells, n_cells)
  /connectivities/{data,indices,indptr}  sparse CSR (n_cells, n_cells)
"""

import sys
from pathlib import Path
from typing import Literal

import anndata as ad
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_knn_parser  # noqa: E402
from schemas import Embedding, Neighbors  # noqa: E402


def compute_neighbors(
    embedding: np.ndarray,
    cell_ids: list[str],
    n_neighbors: int,
    flavor: Literal["umap", "gauss"],
    random_seed: int,
) -> Neighbors:
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


def main() -> None:
    args = build_knn_parser().parse_args()

    emb = Embedding.read(args.pcas)
    nbrs = compute_neighbors(
        emb.matrix,
        emb.row_ids,
        n_neighbors=args.n_neighbors,
        flavor=args.flavor,
        random_seed=args.random_seed,
    )

    out = Path(args.output_dir) / f"{args.name}_neighbors.h5"
    nbrs.write(out)


if __name__ == "__main__":
    main()
