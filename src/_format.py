"""Format-dispatch mixin used by schemas to route ``read``/``write`` to per-format methods."""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

PathLike = str | Path


class FormatDispatch:
    """Mixin: dispatch ``read`` / ``write`` to ``_read_{format}`` / ``_write_{format}``.

    Subclasses must set ``DEFAULT_FORMAT`` and implement ``_write_<fmt>`` /
    ``_read_<fmt>`` for every supported ``fmt``. Unknown formats raise ``ValueError``.
    """

    DEFAULT_FORMAT: ClassVar[str]

    def write(self, path: PathLike, format: str | None = None) -> None:
        fmt = format or self.DEFAULT_FORMAT
        fn = getattr(self, f"_write_{fmt}", None)
        if fn is None:
            raise ValueError(f"unsupported format: {fmt!r}")
        fn(path)

    @classmethod
    def read(cls, path: PathLike, format: str | None = None):
        fmt = format or cls.DEFAULT_FORMAT
        fn = getattr(cls, f"_read_{fmt}", None)
        if fn is None:
            raise ValueError(f"unsupported format: {fmt!r}")
        return fn(path)
