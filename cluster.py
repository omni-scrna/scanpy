#!/usr/bin/env python3
"""Clustering module (scanpy-backed) for omnibenchmark.

Reads UMAP connectivities from HDF5, runs Leiden community detection,
and writes a TSV of cluster assignments.

Output:
  {name}_clusters.tsv  — two-column TSV (cell_id, cluster)
"""

import sys
from pathlib import Path

import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_cluster_parser  # noqa: E402
from readers import read_graph_as_anndata  # noqa: E402
from writers import Clustering, write_clustering  # noqa: E402


def cluster_leiden(adata, resolution=1.0, random_seed=0):
    """Run Leiden (igraph flavor) on ``adata.obsp['connectivities']``."""
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=random_seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata.obs["leiden"]


def main():
    args = build_cluster_parser().parse_args()

    adata, cell_ids = read_graph_as_anndata(args.connectivities)
    labels = cluster_leiden(
        adata, resolution=args.resolution, random_seed=args.random_seed
    )

    out = Path(args.output_dir) / f"{args.name}_clusters.tsv"
    write_clustering(Clustering(cell_ids=cell_ids, labels=labels), out)


if __name__ == "__main__":
    main()
