"""Replay minimum-energy selection for the archived IBM Fez chr22 repeat 1."""
import gzip
import json
from pathlib import Path

import numpy as np


def main():
    with gzip.open(Path(__file__).with_name("ibm_fez_chr22_repeat1.json.gz"), "rt") as f:
        components = json.load(f)
    total = 0.0
    shots = 0
    for component in components:
        counts = component["counts"]
        spins = np.array([[1 if b == "0" else -1 for b in bits[::-1]]
                          for bits in sorted(counts)], dtype=np.int8)
        energies = np.zeros(len(spins))
        for i, j, weight in component["edges"]:
            energies -= weight * (spins[:, i] * spins[:, j])
        total += energies.min()
        shots += sum(counts.values())
    print(f"{len(components)} components, {shots} shots, selected E = {total:.6f}")


if __name__ == "__main__":
    main()
