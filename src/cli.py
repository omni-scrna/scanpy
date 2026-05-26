"""Argument parsers for omnibenchmark scanpy modules.

Reusable across stage entrypoints (pca.py, and future umap.py, cluster.py, ...).

Conventions:
- All arguments are required. No defaults — callers (omnibenchmark configs)
  must pass everything explicitly so runs are fully reproducible from the
  invocation line.
- argparse.parse_args() rejects unknown flags by default; we rely on that
  for strictness rather than parse_known_args().
"""

import argparse


def add_common_args(parser):
    """Args required by omnibenchmark for every module."""
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory for results")
    parser.add_argument("--name", type=str, required=True,
                        help="Module name/identifier")


def build_pca_parser():
    parser = argparse.ArgumentParser(description="OmniBenchmark PCA module (scanpy)")
    add_common_args(parser)

    parser.add_argument("--normalized_selected.h5", dest="input_h5",
                        type=str, required=True,
                        help="TENx-format HDF5 of normalized, selected expression (genes x cells)")

    parser.add_argument("--solver", type=str, required=False, default=None,
                        choices=["arpack", "randomized"],
                        help="PCA solver. Required when --chunked=false; must be omitted when --chunked=true (IncrementalPCA ignores it).")
    parser.add_argument("--n_components", type=int, required=True,
                        help="Number of principal components to compute")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Seed for randomized solvers (and for reproducibility)")
    parser.add_argument("--chunked", type=str, required=True,
                        choices=["true", "false"],
                        help="If true, use incremental PCA (IncrementalPCA) over chunks")
    parser.add_argument("--chunk_size", type=int, required=False, default=None,
                        help="Chunk size for incremental PCA. Required when --chunked=true; must be omitted otherwise.")

    return parser


def build_knn_parser():
    parser = argparse.ArgumentParser(description="OmniBenchmark kNN module (scanpy)")
    add_common_args(parser)

    parser.add_argument("--pca.tsv", dest="pca_tsv", type=str, required=True,
                        help="PCA embedding TSV produced by the pca entrypoint")
    parser.add_argument("--n_neighbors", type=int, required=True,
                        help="Number of nearest neighbors")
    parser.add_argument("--flavor", type=str, required=True,
                        choices=["umap", "gauss"],
                        help="Method to compute connectivities")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Random seed")

    return parser
