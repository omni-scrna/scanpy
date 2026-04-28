#!/usr/bin/env python3
# test that all runtime dependencies import cleanly
import h5py          # noqa: F401
import numpy         # noqa: F401
import scipy.sparse  # noqa: F401
import anndata       # noqa: F401
import scanpy        # noqa: F401

print("OK")
