#!/usr/bin/env python3
"""kNN graph module (scanpy-backed) for omnibenchmark.

Output HDF5: {output_dir}/{name}_neighbors.h5
  Flat CSR of the kNN distance graph, matching the R metrics reader
  (graph.R::read_csr_h5), which reads these datasets from the file root:
    /cell_ids   string array (n_cells,)
    /data       CSR data
    /indices    CSR column indices (0-based)
    /indptr     CSR row pointers
"""

import argparse
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import polars as pl
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))  # vendored `common` package (src/common)
from common import cli  # noqa: E402

def parse_args():
    # We own the parser; common/cli injects the shared contract (base args + the
    # `nngraph` stage I/O from common/schema). This module's method params are
    # hand-rolled below, so the whole CLI stays visible here.
    p = argparse.ArgumentParser(description="kNN graph module (scanpy-backed)")
    cli.add_base_args(p)                # --output_dir, --name
    cli.add_stage_args(p, "NNG")    # --pcas_tsv
    p.add_argument("--n_neighbors", type=int, required=True,
                   help="Number of nearest neighbors")
    p.add_argument("--flavor", type=str, required=True,
                   choices=["umap", "gauss"], help="Method to compute connectivities")
    p.add_argument("--random_seed", type=int, required=True, help="Random seed")
    return p.parse_args()


def _write_sparse(h5, cell_ids, m):
    # Flat layout with datasets at the file root.
    m = m.tocsr()
    # dtype="S": h5py can't write numpy unicode ('<U') arrays; bytes give portable fixed-length HDF5 strings.
    h5.create_dataset("cell_ids", data=np.array(cell_ids, dtype="S"))
    h5.create_dataset("data",    data=m.data)
    h5.create_dataset("indices", data=m.indices)
    h5.create_dataset("indptr",  data=m.indptr)


def write_neighbors_graph(adata, out_dir, name):
    out = Path(out_dir) / f"{name}_neighbors.h5"
    with h5py.File(out, "w") as h5:
        _write_sparse(h5, adata.obs_names.to_list(), adata.obsp["distances"])
    print(f"  wrote: {out}")


def main():
    args = parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "pcas_tsv", "n_neighbors", "flavor", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # TSV has N header cols and N+1 data cols (first data col = row IDs, unnamed).
    df = pl.read_csv(args.pcas_tsv, separator="\t", skip_rows=1, has_header=False)
    embedding = df[:, 1:].to_numpy().astype(np.float64)

    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
    adata.obs_names = df[:, 0].to_list()
    adata.obsm["X_pca"] = embedding

    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, method=args.flavor,
                    use_rep="X_pca", random_state=args.random_seed)

    write_neighbors_graph(adata, args.output_dir, args.name)


if __name__ == "__main__":
    main()
