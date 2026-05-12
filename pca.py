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

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent / "src"))
from cli import build_pca_parser  # noqa: E402
from readers import read_tenx_h5  # noqa: E402
from schemas import Embedding  # noqa: E402


def run_pca(adata: ad.AnnData, args: argparse.Namespace) -> Embedding:
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
    # the data is sparse, Scanpy may need to densify each chunk during the
    # partial_fit step, which can cause local spikes in memory usage.

    sc.pp.pca(
        adata,
        n_comps=args.n_components,
        zero_center=True,
        svd_solver=args.solver,
        random_state=args.random_seed,
    )

    matrix = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    return Embedding(
        matrix=matrix,
        row_ids=list(adata.obs_names),
        col_names=[f"PC{i + 1}" for i in range(matrix.shape[1])],
    )


def main() -> None:
    args = build_pca_parser().parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    adata = read_tenx_h5(args.input_h5)
    print(f"  matrix (cells x genes): {adata.shape}")

    emb = run_pca(adata, args)
    print(f"  embedding: {emb.matrix.shape}")

    out = Path(args.output_dir) / f"{args.name}_pcas.tsv"
    emb.write(out)
    print(f"  wrote: {out}")


if __name__ == "__main__":
    main()
