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
- Genes are always centered/scaled before PCA (sc.pp.scale, zero_center=True).
  If alternative scaling is needed later, maybe expose it as a new --pca_type
  variant rather than as an independent flag.
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
    sc.pp.scale(adata, zero_center=True, max_value=None)

    # TODO: accept chunked argument too, to cap memory usage (triggers
    # incremental PCA approach -- maybe pass it as pca_type and enforce
    # chunked.
    # A larger chunk_size improves numerical stability and speed (up to a
    # point) but increases the memory ceiling.
    # Chunked mode is most effective when paired with backed mode (e.g., adata
    # = sc.read_h5ad('file.h5ad', backed='r')). If the AnnData object is
    # already fully loaded into memory, using chunked=True provides no memory
    # benefit and will simply slow down the computation.
    # IncrementalPCA in scikit-learn traditionally prefers dense inputs. If
    # your data is sparse, Scanpy may need to densify each chunk during the
    # partial_fit step, which can cause local spikes in memory usage.

    sc.pp.pca(
        adata,
        n_comps=args.n_components,
        zero_center=True,
        svd_solver=args.solver,
        random_state=args.random_seed,
    )

    embedding = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    loadings = np.asarray(adata.varm["PCs"], dtype=np.float64)
    variance = np.asarray(adata.uns["pca"]["variance"], dtype=np.float64)
    variance_ratio = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=np.float64)
    return embedding, loadings, variance, variance_ratio



def main():
    args = build_pca_parser().parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "input_h5", "solver", "n_components", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    adata = load_matrix(args.input_h5)
    gene_ids = np.array(adata.var_names)
    cell_ids = np.array(adata.obs_names)
    print(f"  matrix (cells x genes): {adata.shape}")

    embedding, loadings, variance, variance_ratio = run_pca(adata, args)
    print(f"  embedding: {embedding.shape}, loadings: {loadings.shape}")

    col_names = [f"PC{i + 1}" for i in range(embedding.shape[1])]
    out = Path(args.output_dir) / f"{args.name}_pcas.tsv"
    write_embeddings(Embedding(embedding, list(cell_ids), col_names), out)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
