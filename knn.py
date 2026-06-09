#!/usr/bin/env python3
"""kNN graph module (scanpy-backed) for omnibenchmark.

Output HDF5: {output_dir}/{name}_knn.h5
  /distances/data, /distances/indices, /distances/indptr, /distances/shape
  /connectivities/data, /connectivities/indices, /connectivities/indptr, /connectivities/shape
  attrs: format_version, tool, tool_version, n_neighbors, random_seed
"""

import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import polars as pl
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))  # vendored `common` package
from common.cli import parse_args  # noqa: E402

OUTPUT_FORMAT_VERSION = "1"
TOOL = "scanpy"


def write_sparse(h5, name, m):
    m = m.tocsr()
    g = h5.create_group(name)
    g.create_dataset("data",    data=m.data)
    g.create_dataset("indices", data=m.indices)
    g.create_dataset("indptr",  data=m.indptr)
    g.create_dataset("shape",   data=np.array(m.shape))


def main():
    args = parse_args("knn")
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "pca_tsv", "n_neighbors", "flavor", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # TSV has N header cols and N+1 data cols (first data col = row IDs, unnamed).
    df = pl.read_csv(args.pca_tsv, separator="\t", skip_rows=1, has_header=False)
    embedding = df[:, 1:].to_numpy().astype(np.float64)

    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
    adata.obs_names = df[:, 0].to_list()
    adata.obsm["X_pca"] = embedding

    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, method=args.flavor,
                    use_rep="X_pca", random_state=args.random_seed)

    out = Path(args.output_dir) / f"{args.name}_knn.h5"
    with h5py.File(out, "w") as h5:
        write_sparse(h5, "distances",     adata.obsp["distances"])
        write_sparse(h5, "connectivities", adata.obsp["connectivities"])
        from importlib.metadata import version
        h5.attrs["format_version"] = OUTPUT_FORMAT_VERSION
        h5.attrs["tool"]           = TOOL
        h5.attrs["tool_version"]   = version("scanpy")
        h5.attrs["n_neighbors"]    = args.n_neighbors
        h5.attrs["flavor"]         = args.flavor
        h5.attrs["random_seed"]    = args.random_seed

    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
