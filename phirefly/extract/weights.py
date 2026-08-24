"""Observation weighting helpers."""

from __future__ import annotations

import math


def observation_weight(
    baseq: int,
    mapq: int,
    epsilon_floor: float,
    epsilon_ceiling: float,
    mapq_cap: int,
    cache: dict[tuple[int, int], tuple[float, float]],
) -> tuple[float, float]:
    key = (int(baseq), min(int(mapq), int(mapq_cap)))
    cached = cache.get(key)
    if cached is not None:
        return cached
    base_error = 10 ** (-key[0] / 10)
    map_error = 10 ** (-key[1] / 10)
    epsilon = 1 - (1 - base_error) * (1 - map_error)
    epsilon = min(max(epsilon, epsilon_floor), epsilon_ceiling)
    weight = math.log((1 - epsilon) / epsilon)
    cache[key] = (epsilon, weight)
    return epsilon, weight


__all__ = ["observation_weight"]
