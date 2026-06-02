#!/usr/bin/env python3
"""
PCA module (scanpy-backed) for omnibenchmark.

Output
------
File: {output_dir}/{name}_pcas.tsv

Tab-separated, with header row:
  cell_id  PC1  PC2  ...  PC{n_components}

One row per cell; values are float64 PCA scores.

Implementation notes
--------------------
- Genes are mean-centered (only) before PCA via sc.pp.pca(zero_center=True).
  No per-gene variance scaling — matches scrapper/rapids-singlecell. If
  alternative scaling is needed later, expose it as a new --pca_type variant
  rather than as an independent flag.
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_pca_parser  # noqa: E402
from writers import Embedding, write_embeddings  # noqa: E402
from phases import phase  # noqa: E402
from obkit.logger import init_logger  # noqa: E402


def load_matrix(h5_path):
    with h5py.File(h5_path, "r") as h5:
        g = h5["matrix"]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])
        gene_ids = g["genes"][:].astype(str)
        cell_ids = g["barcodes"][:].astype(str)

    X = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()  # cells x genes
    adata = ad.AnnData(X=X)
    adata.obs_names = cell_ids
    adata.var_names = gene_ids
    return adata


def run_pca(adata, args):
    # Chunked mode triggers IncrementalPCA. It is most effective with backed
    # AnnData; for in-memory data it provides no memory benefit and is slower.
    # Sparse inputs get densified per chunk during partial_fit.
    chunked = args.chunked == "true"

    sc.pp.pca(
        adata,
        n_comps=args.n_components,
        zero_center=True,
        svd_solver=args.solver if not chunked else None,
        random_state=args.random_seed,
        chunked=chunked,
        chunk_size=args.chunk_size if chunked else None,
    )

    embedding = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    loadings = np.asarray(adata.varm["PCs"], dtype=np.float64)
    variance = np.asarray(adata.uns["pca"]["variance"], dtype=np.float64)
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=np.float64)
    return embedding, loadings, variance, variance_ratio


def validate_args(args):
    chunked = args.chunked == "true"
    if chunked:
        if args.chunk_size is None:
            sys.exit("error: --chunk_size is required when --chunked=true")
        if args.solver is not None:
            sys.exit("error: --solver must not be set when --chunked=true (IncrementalPCA ignores it)")
    else:
        if args.chunk_size is not None:
            sys.exit("error: --chunk_size must not be set when --chunked=false")
        if args.solver is None:
            sys.exit("error: --solver is required when --chunked=false")



def main():
    args = build_pca_parser().parse_args()
    validate_args(args)
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "input_h5", "solver", "n_components", "random_seed", "chunked", "chunk_size"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    init_logger(args.output_dir)

    with phase("load") as attrs:
        adata = load_matrix(args.input_h5)
        attrs["n_cells"] = adata.n_obs
        attrs["n_genes"] = adata.n_vars
    gene_ids = np.array(adata.var_names)
    cell_ids = np.array(adata.obs_names)
    print(f"  matrix (cells x genes): {adata.shape}")

    with phase("pca") as attrs:
        embedding, loadings, variance, variance_ratio = run_pca(adata, args)
        attrs["solver"] = args.solver or "chunked"
        attrs["n_components"] = args.n_components
    print(f"  embedding: {embedding.shape}, loadings: {loadings.shape}")

    with phase("write"):
        col_names = [f"PC{i + 1}" for i in range(embedding.shape[1])]
        out = Path(args.output_dir) / f"{args.name}_pcas.tsv"
        write_embeddings(Embedding(embedding, list(cell_ids), col_names), out)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
