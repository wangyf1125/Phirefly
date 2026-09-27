"""Small file I/O helpers shared across Phirefly modules."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    path_s = str(path)
    if path_s.endswith(".gz"):
        return gzip.open(path_s, mode)
    return open(path_s, mode)


def read_key_value(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with open_text(path) as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        return {row[0]: row[1] for row in rows if len(row) >= 2}


def require_fresh_output(path: str | Path, force: bool = False) -> None:
    """Never reuse partial or stale results implicitly."""
    path = Path(path)
    if path.exists() and not path.is_dir():
        raise ValueError(f"Output path is not a directory: {path}")
    if path.is_dir() and any(path.iterdir()) and not force:
        raise ValueError(f"Output directory is not empty: {path}; use a new directory or --force to recompute")


__all__ = ["open_text", "read_key_value", "require_fresh_output"]
