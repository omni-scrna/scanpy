#!/usr/bin/env python3
"""
PCA module (scanpy-backed) for omnibenchmark.

Output format (neutral HDF5, version 1)
---------------------------------------
File: {output_dir}/{name}_pca.h5

Datasets:
  /embedding       float64, shape (n_cells, n_components)
                   PCA scores. Cells-as-rows. NOTE for R consumers:
                   transpose with t() if you need cells-as-columns.
                   (SingleCellExperiment reducedDim is also cells-as-rows,
                   so usually no transpose is needed.)
  /loadings        float64, shape (n_genes, n_components)
                   Gene loadings (a.k.a. rotation / components). Genes-as-rows.
  /variance        float64, shape (n_components,)
                   Variance explained by each PC.
  /variance_ratio  float64, shape (n_components,)
                   Fraction of total variance explained per PC.
  /cell_ids        utf-8 string, shape (n_cells,)
                   Cell barcodes, aligned with rows of /embedding.
  /gene_ids        utf-8 string, shape (n_genes,)
                   Gene identifiers (the selected/HVG subset),
                   aligned with rows of /loadings.

Root attributes (provenance / params):
  tool="scanpy", tool_version, solver,
  n_components, random_seed, format_version="1"

This format is intentionally flat HDF5 (not h5ad) so R (rhdf5 / HDF5Array)
and Julia (HDF5.jl) implementations can read/write it without anndata.

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

OUTPUT_FORMAT_VERSION = "1"
TOOL = "scanpy"


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


def write_output(path, embedding, loadings, variance, variance_ratio,
                 cell_ids, gene_ids, args):
    from importlib.metadata import version

    str_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(path, "w") as h5:
        h5.create_dataset("embedding", data=embedding)
        h5.create_dataset("loadings", data=loadings)
        h5.create_dataset("variance", data=variance)
        h5.create_dataset("variance_ratio", data=variance_ratio)
        h5.create_dataset("cell_ids", data=np.asarray(cell_ids, dtype=object), dtype=str_dtype)
        h5.create_dataset("gene_ids", data=np.asarray(gene_ids, dtype=object), dtype=str_dtype)

        h5.attrs["format_version"] = OUTPUT_FORMAT_VERSION
        h5.attrs["tool"] = TOOL
        h5.attrs["tool_version"] = version("scanpy")
        h5.attrs["solver"] = args.solver
        h5.attrs["n_components"] = args.n_components
        h5.attrs["random_seed"] = args.random_seed


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

    out = Path(args.output_dir) / f"{args.name}_pca.h5"
    write_output(out, embedding, loadings, variance, variance_ratio,
                 cell_ids, gene_ids, args)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
