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

import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_pca_parser  # noqa: E402
from readers import read_tenx_h5  # noqa: E402
from schemas import Embedding  # noqa: E402


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

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    adata = read_tenx_h5(args.input_h5)
    cell_ids = np.array(adata.obs_names)
    print(f"  matrix (cells x genes): {adata.shape}")

    embedding, loadings, variance, variance_ratio = run_pca(adata, args)
    print(f"  embedding: {embedding.shape}, loadings: {loadings.shape}")

    col_names = [f"PC{i + 1}" for i in range(embedding.shape[1])]
    out = Path(args.output_dir) / f"{args.name}_pcas.tsv"
    Embedding(embedding, list(cell_ids), col_names).write(out)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
