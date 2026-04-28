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

    parser.add_argument("--normalized.h5", dest="normalized_h5",
                        type=str, required=True,
                        help="TENx-format HDF5 of normalized expression (genes x cells)")
    parser.add_argument("--selected.genes", dest="selected_genes",
                        type=str, required=True,
                        help="Gzipped text file of selected gene ids (one per line)")

    parser.add_argument("--pca_type", type=str, required=True,
                        choices=["scanpy_arpack", "scanpy_randomized"],
                        help="PCA solver / implementation")
    parser.add_argument("--n_components", type=int, required=True,
                        help="Number of principal components to compute")
    parser.add_argument("--random_seed", type=int, required=True,
                        help="Seed for randomized solvers (and for reproducibility)")

    return parser
