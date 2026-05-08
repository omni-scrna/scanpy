#!/usr/bin/env python3
"""Clustering module (scanpy-backed) for omnibenchmark.

Reads UMAP connectivities from HDF5, runs Leiden community detection,
and writes a TSV of cluster assignments.

Output:
  {name}_clusters.tsv  — two-column TSV (cell_id, cluster)
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_cluster_parser  # noqa: E402
from writers import Clustering, read_graph, write_clustering  # noqa: E402


def main():
    args = build_cluster_parser().parse_args()

    graph = read_graph(args.connectivities)
    n_cells = len(graph.cell_ids)

    adata = ad.AnnData(X=np.zeros((n_cells, 1)))
    adata.obs_names = graph.cell_ids
    adata.obsp["connectivities"] = graph.matrix
    adata.uns["neighbors"] = {"connectivities_key": "connectivities"}

    sc.tl.leiden(adata, resolution=args.resolution, random_state=args.random_seed)

    out = Path(args.output_dir) / f"{args.name}_clusters.tsv"
    write_clustering(
        Clustering(cell_ids=graph.cell_ids, labels=adata.obs[cluster_col].values), out
    )


if __name__ == "__main__c":
    main()
