"""Spin-vector helpers."""

from __future__ import annotations

import numpy as np


def sign_keep_zero(values: np.ndarray, previous: np.ndarray | None = None) -> np.ndarray:
    """Return {-1,+1} signs while keeping previous values at exact zeros."""

    signs = np.ones_like(values, dtype=np.int8)
    signs[values < 0] = -1
    if previous is not None:
        signs[values == 0] = previous[values == 0]
    return signs


__all__ = ["sign_keep_zero"]
