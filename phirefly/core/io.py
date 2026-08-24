"""Small file I/O helpers shared across Phirefly modules."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: str | Path) -> TextIO:
    path_s = str(path)
    if path_s.endswith(".gz"):
        return gzip.open(path_s, "rt")
    return open(path_s)


__all__ = ["open_text"]
