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
import polars as pl
import scipy.sparse as sp
import yaml
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



def read_batch_assignments(rawdata_h5ad, batch_variable):
    # Read batch labels from h5ad obs
    with h5py.File(rawdata_h5ad, "r") as f:
        cell_ids = f["obs/_index"][:].astype(str)
        batch_data = f[f"obs/{batch_variable}"]
        if "categories" in batch_data:
            codes = batch_data["codes"][:]
            categories = batch_data["categories"][:].astype(str)
            batch = categories[codes]
        else:
            batch = batch_data[:].astype(str)
    return cell_ids, batch


def run_pca_per_batch(adata, args):
    # compute PCA per batch, return concatenated embeddings and per-batch loadings
    with open(args.batch_info_yaml) as f:
        batch_info = yaml.safe_load(f)
    batch_variable = batch_info["batch_var"]

    cell_ids_all, batch_all = read_batch_assignments(args.rawdata_h5ad, batch_variable)
    cell_ids_matrix = np.array(adata.obs_names)

    idx = np.array([np.where(cell_ids_all == c)[0][0] for c in cell_ids_matrix])
    batch_aligned = batch_all[idx]
    batches = np.unique(batch_aligned)
    print(f"  per-batch PCA: {len(batches)} batches ({', '.join(batches)})")

    gene_names = np.array(adata.var_names)
    embeddings_parts = []
    loadings_dict = {}

    for b in batches:
        mask = batch_aligned == b
        adata_batch = adata[mask].copy()
        print(f"    batch '{b}': {adata_batch.n_obs} cells")

        embedding, loadings, _, _ = run_pca(adata_batch, args)
        col_names = [f"PC{i + 1}" for i in range(embedding.shape[1])]

        batch_df = pl.DataFrame({
            "cell_id": np.array(adata_batch.obs_names).tolist(),
            **{col: embedding[:, i].tolist() for i, col in enumerate(col_names)},
            "batch_id": [b] * adata_batch.n_obs,
        })
        embeddings_parts.append(batch_df)
        loadings_dict[b] = loadings

    return pl.concat(embeddings_parts), loadings_dict, gene_names


def main():
    args = build_pca_parser().parse_args()
    print(f"Full command: {' '.join(sys.argv)}")
    for k in ("output_dir", "name", "input_h5", "solver", "n_components", "random_seed", "per_batch"):
        print(f"  {k}: {getattr(args, k)}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    adata = load_matrix(args.input_h5)
    gene_ids = np.array(adata.var_names)
    cell_ids = np.array(adata.obs_names)
    print(f"  matrix (cells x genes): {adata.shape}")

    if args.per_batch == "true":
        embeddings_df, loadings_dict, gene_names = run_pca_per_batch(adata, args)

        out_tsv = Path(args.output_dir) / f"{args.name}_pcas_per_batch.tsv"
        embeddings_df.write_csv(out_tsv, separator="\t")
        print(f"  wrote: {out_tsv}")

        out_h5 = Path(args.output_dir) / f"{args.name}_loadings_per_batch.h5"
        with h5py.File(out_h5, "w") as f:
            f.create_dataset("gene_names", data=gene_names.astype(bytes))
            for b, loadings in loadings_dict.items():
                grp = f.create_group(b)
                grp.create_dataset("loadings", data=loadings)
        print(f"  wrote: {out_h5}")
    else:
        embedding, loadings, variance, variance_ratio = run_pca(adata, args)
        print(f"  embedding: {embedding.shape}, loadings: {loadings.shape}")

        col_names = [f"PC{i + 1}" for i in range(embedding.shape[1])]
        out = Path(args.output_dir) / f"{args.name}_pcas.tsv"
        write_embeddings(Embedding(embedding, list(cell_ids), col_names), out)
        print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
