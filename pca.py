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
  pca_type, n_components, random_seed,
  scanpy_version, anndata_version, format_version="1"

This format is intentionally flat HDF5 (not h5ad) so R (rhdf5 / HDF5Array)
and Julia (HDF5.jl) implementations can read/write it without anndata.

Implementation notes
--------------------
- Genes are always centered/scaled before PCA (sc.pp.scale, zero_center=True).
  If alternative scaling is needed later, maybe expose it as a new --pca_type
  variant rather than as an independent flag.
- Subsetting to selected genes happens here for now; this responsibility
  should move to a dedicated upstream cleanup stage. See load_subset_matrix.
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_pca_parser  # noqa: E402

OUTPUT_FORMAT_VERSION = "1"

SOLVER = {"scanpy_arpack": "arpack", "scanpy_randomized": "randomized"}


def _decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def load_selected_genes(path):
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def load_subset_matrix(h5_path, selected_genes):
    """Load TENx HDF5 (genes x cells) and subset to selected_genes.

    Returns (X_csr_cells_x_genes, gene_ids, cell_ids).
    """
    with h5py.File(h5_path, "r") as h5:
        g = h5["matrix"]
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])  # (n_genes, n_cells)
        if "features" in g and "id" in g["features"]:
            gene_ids = g["features/id"][:]
        elif "genes" in g:
            gene_ids = g["genes"][:]
        else:
            gene_ids = np.array([f"gene_{i}".encode() for i in range(shape[0])])
        cell_ids = g["barcodes"][:] if "barcodes" in g else \
            np.array([f"cell_{i}".encode() for i in range(shape[1])])

    gene_ids = _decode(gene_ids)
    cell_ids = _decode(cell_ids)

    # TENx HDF5 stores data in CSC order (indptr = column pointers).
    m = sp.csc_matrix((data, indices, indptr), shape=shape)

    sel_set = set(selected_genes)
    available = set(gene_ids.tolist())
    missing = sel_set - available
    if missing:
        sample = sorted(missing)[:10]
        raise ValueError(
            f"{len(missing)} selected gene(s) not present in normalized.h5; "
            f"first {len(sample)}: {sample}"
        )
    mask = np.array([gene in sel_set for gene in gene_ids])

    m_sub = m[mask, :]
    X = m_sub.T.tocsr()  # cells-as-rows
    return X, gene_ids[mask], cell_ids


def run_pca(X, gene_ids, cell_ids, args):
    import anndata as ad
    import scanpy as sc

    adata = ad.AnnData(X=X)
    adata.obs_names = cell_ids
    adata.var_names = gene_ids

    sc.pp.scale(adata, zero_center=True, max_value=None)

    solver = SOLVER[args.pca_type]

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
        svd_solver=solver,
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
        h5.attrs["pca_type"] = args.pca_type
        h5.attrs["n_components"] = args.n_components
        h5.attrs["random_seed"] = args.random_seed
        h5.attrs["scanpy_version"] = version("scanpy")
        h5.attrs["anndata_version"] = version("anndata")


def main():
    args = build_pca_parser().parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "normalized_h5", "selected_genes",
              "pca_type", "n_components", "random_seed"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    selected = load_selected_genes(args.selected_genes)
    print(f"  selected genes: {len(selected)}")

    X, gene_ids, cell_ids = load_subset_matrix(args.normalized_h5, selected)
    print(f"  matrix (cells x genes): {X.shape}")

    embedding, loadings, variance, variance_ratio = run_pca(X, gene_ids, cell_ids, args)
    print(f"  embedding: {embedding.shape}, loadings: {loadings.shape}")

    out = Path(args.output_dir) / f"{args.name}_pca.h5"
    write_output(out, embedding, loadings, variance, variance_ratio,
                 cell_ids, gene_ids, args)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
