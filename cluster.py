#!/usr/bin/env python3
"""Leiden clustering module (scanpy-backed) for omnibenchmark.

Reads {name}_neighbors.h5 (the kNN distances + connectivities from the knn
entrypoint), runs Leiden on the connectivities, and writes:
  {output_dir}/{name}_clusters.tsv  — cell_id<TAB>cluster
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import polars as pl
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))  # vendored `common` package (src/common)
from common import cli  # noqa: E402
from readers import read_neighbors  # noqa: E402


def parse_args():
    # common/cli injects the shared contract (base args + the CLUST stage I/O);
    # method params are hand-rolled below.
    p = argparse.ArgumentParser(description="Leiden clustering module (scanpy-backed)")
    cli.add_base_args(p)             # --output_dir, --name
    cli.add_stage_args(p, "CLUST")   # --neighbors_h5
    p.add_argument("--resolution", type=float, required=True,
                   help="Resolution controlling cluster granularity")
    p.add_argument("--random_seed", type=int, required=True, help="Random seed")
    return p.parse_args()


def build_adata(neighbors_h5):
    """AnnData with the stored distances/connectivities graph wired up for scanpy."""
    distances, connectivities, cell_ids = read_neighbors(neighbors_h5)
    adata = ad.AnnData(X=np.zeros((len(cell_ids), 1)))
    adata.obs_names = cell_ids
    adata.obsp["distances"] = distances
    adata.obsp["connectivities"] = connectivities
    adata.uns["neighbors"] = {
        "distances_key": "distances",
        "connectivities_key": "connectivities",
    }
    return adata, cell_ids


def cluster_leiden(adata, resolution, random_seed):
    sc.tl.leiden(adata, resolution=resolution, random_state=random_seed,
                 flavor="igraph", n_iterations=2, directed=False)
    return adata.obs["leiden"].to_list()


def main():
    args = parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "neighbors_h5", "resolution", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    adata, cell_ids = build_adata(args.neighbors_h5)
    labels = cluster_leiden(adata, args.resolution, args.random_seed)

    out = Path(args.output_dir) / f"{args.name}_clusters.tsv"
    pl.DataFrame({"cell_id": cell_ids, "cluster": labels}).write_csv(out, separator="\t")
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
