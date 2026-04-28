#!/usr/bin/env python3
"""
Validate a PCA output HDF5 file against the format-version 1 spec.

Usage
-----
    python validators/pca_output.py path/to/name_pca.h5

Exit codes: 0 = valid, 1 = validation failure, 2 = file not found / IO error.
"""

import sys
from pathlib import Path

import h5py
import numpy as np


REQUIRED_DATASETS = {
    "embedding":       ("float64", 2),
    "loadings":        ("float64", 2),
    "variance":        ("float64", 1),
    "variance_ratio":  ("float64", 1),
    "cell_ids":        (None,      1),
    "gene_ids":        (None,      1),
}

REQUIRED_ATTRS = {"format_version", "pca_type", "n_components", "random_seed",
                  "scanpy_version", "anndata_version"}

VALID_PCA_TYPES = {"scanpy_arpack", "scanpy_randomized"}


def validate(path) -> list[str]:
    """Return a list of error strings; empty list means the file is valid."""
    errors = []
    path = Path(path)

    if not path.exists():
        return [f"file not found: {path}"]

    try:
        h5 = h5py.File(path, "r")
    except OSError as e:
        return [f"cannot open HDF5 file: {e}"]

    with h5:
        # ── root attributes ──────────────────────────────────────────────────
        missing_attrs = REQUIRED_ATTRS - set(h5.attrs)
        for a in sorted(missing_attrs):
            errors.append(f"missing root attribute: {a}")

        if "format_version" in h5.attrs and h5.attrs["format_version"] != "1":
            errors.append(
                f"unexpected format_version: {h5.attrs['format_version']!r} (expected '1')"
            )

        if "pca_type" in h5.attrs and h5.attrs["pca_type"] not in VALID_PCA_TYPES:
            errors.append(
                f"unknown pca_type: {h5.attrs['pca_type']!r} "
                f"(valid: {sorted(VALID_PCA_TYPES)})"
            )

        # ── datasets: existence, dtype, ndim ────────────────────────────────
        for name, (expected_dtype, expected_ndim) in REQUIRED_DATASETS.items():
            if name not in h5:
                errors.append(f"missing dataset: /{name}")
                continue
            ds = h5[name]
            if expected_dtype and ds.dtype != np.dtype(expected_dtype):
                errors.append(
                    f"/{name}: expected dtype {expected_dtype}, got {ds.dtype}"
                )
            if ds.ndim != expected_ndim:
                errors.append(
                    f"/{name}: expected {expected_ndim}D array, got {ds.ndim}D"
                )

        # ── cross-dataset shape consistency ──────────────────────────────────
        if all(k in h5 for k in ("embedding", "loadings", "variance",
                                  "variance_ratio", "cell_ids", "gene_ids")):
            n_cells  = h5["embedding"].shape[0]
            n_comps  = h5["embedding"].shape[1]
            n_genes  = h5["loadings"].shape[0]

            if h5["loadings"].shape[1] != n_comps:
                errors.append(
                    f"/loadings second dim {h5['loadings'].shape[1]} != "
                    f"n_components {n_comps} from /embedding"
                )
            if h5["variance"].shape[0] != n_comps:
                errors.append(
                    f"/variance length {h5['variance'].shape[0]} != n_components {n_comps}"
                )
            if h5["variance_ratio"].shape[0] != n_comps:
                errors.append(
                    f"/variance_ratio length {h5['variance_ratio'].shape[0]} != "
                    f"n_components {n_comps}"
                )
            if h5["cell_ids"].shape[0] != n_cells:
                errors.append(
                    f"/cell_ids length {h5['cell_ids'].shape[0]} != n_cells {n_cells}"
                )
            if h5["gene_ids"].shape[0] != n_genes:
                errors.append(
                    f"/gene_ids length {h5['gene_ids'].shape[0]} != n_genes {n_genes}"
                )

            if "n_components" in h5.attrs:
                stored_n_comps = int(h5.attrs["n_components"])
                if stored_n_comps != n_comps:
                    errors.append(
                        f"n_components attribute ({stored_n_comps}) != "
                        f"actual embedding columns ({n_comps})"
                    )

        # ── value sanity checks ───────────────────────────────────────────────
        if "variance_ratio" in h5:
            vr = h5["variance_ratio"][:]
            if np.any(vr < 0) or np.any(vr > 1):
                errors.append("/variance_ratio values must be in [0, 1]")
            if vr.sum() > 1.0 + 1e-6:
                errors.append(
                    f"/variance_ratio sums to {vr.sum():.6f}, expected <= 1.0"
                )

        if "variance" in h5:
            v = h5["variance"][:]
            if np.any(v < 0):
                errors.append("/variance contains negative values")

    return errors


def main():
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <pca_output.h5>", file=sys.stderr)
        sys.exit(2)

    path = sys.argv[1]
    errors = validate(path)

    if errors:
        print(f"INVALID: {path}")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print(f"OK: {path}")
        with h5py.File(path, "r") as h5:
            n_cells = h5["embedding"].shape[0]
            n_comps = h5["embedding"].shape[1]
            n_genes = h5["loadings"].shape[0]
            print(f"  cells={n_cells}  genes={n_genes}  components={n_comps}")
            print(f"  pca_type={h5.attrs['pca_type']}  "
                  f"scanpy={h5.attrs['scanpy_version']}")
        sys.exit(0)


if __name__ == "__main__":
    main()
