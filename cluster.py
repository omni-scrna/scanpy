#!/usr/bin/env python3
"""Clustering module (scanpy-backed) for omnibenchmark.

Reads a Neighbors bundle from HDF5 (see ``knn.py`` for the layout),
runs Leiden community detection on the connectivities, and writes a
TSV of cluster assignments.

Output:
  {output_dir}/{name}_clusters.tsv — two-column TSV (cell_id, cluster)
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_cluster_parser  # noqa: E402
from readers import read_neighbors_as_anndata  # noqa: E402
from schemas import Clustering  # noqa: E402


def cluster_leiden(
    adata: ad.AnnData, resolution: float = 1.0, random_seed: int = 0
) -> pd.Series:
    """Run Leiden (igraph flavor) on ``adata.obsp['connectivities']``."""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return pd.Series(adata.obs["leiden"])


def main() -> None:
    args = build_cluster_parser().parse_args()

    adata, cell_ids = read_neighbors_as_anndata(args.neighbors)
    labels = cluster_leiden(
        adata, resolution=args.resolution, random_seed=args.random_seed
    )

    out = Path(args.output_dir) / f"{args.name}_clusters.tsv"
    Clustering(cell_ids=cell_ids, labels=list(labels)).write(out)


if __name__ == "__main__":
    main()
