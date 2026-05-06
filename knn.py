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
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_knn_parser  # noqa: E402

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
    args = build_knn_parser().parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "pca_tsv", "n_neighbors", "flavor", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    df = pd.read_table(args.pca_tsv, index_col=0)
    embedding = df.to_numpy(dtype=np.float64)

    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 1)))
    adata.obs_names = df.index.tolist()
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
