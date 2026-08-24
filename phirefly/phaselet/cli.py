#!/usr/bin/env python3
"""CLI compatibility wrapper for the phaselet/hyperread Phirefly-QAIA runner."""

from __future__ import annotations

import sys
from pathlib import Path

try:
    from .runner import main
except ImportError:  # pragma: no cover - direct script execution fallback
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from phirefly.phaselet.runner import main  # type: ignore


if __name__ == "__main__":
    main()
