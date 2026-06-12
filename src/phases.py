"""Phase boundary helper around obkit.logger.

Copied from the rapids-singlecell module to keep wire format identical.
TODO: promote to obkit so both modules can `from obkit.logger import phase`.

Usage:

    init_logger(output_dir)
    with phase("load") as attrs:
        adata = load(...)
        attrs["n_cells"] = adata.n_obs
"""

from contextlib import contextmanager

from obkit.logger import emit


@contextmanager
def phase(name):
    attrs = {}
    emit(name, "start")
    try:
        yield attrs
    finally:
        emit(name, "end", attrs=attrs or None)
