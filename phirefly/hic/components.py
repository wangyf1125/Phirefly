"""Component and SNP-site helpers for Hi-C/Pore-C orientation edges."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from ..core.genomics import parse_snp_id
from ..vcf import load_component_phase_sets, load_phases


@dataclass(frozen=True)
class ComponentSite:
    ref: str
    alt: str
    delta: int
    component: str
    snp_id: str


def load_component_site_map(
    observations: str | Path,
    snp_phases: str | Path,
) -> tuple[dict[str, dict[int, list[ComponentSite]]], dict[str, list[int]]]:
    """Build chrom/position lookup for phased SNPs and their read-SNP components."""

    phases = load_phases(snp_phases)
    components = load_component_phase_sets(observations)
    site_map: dict[str, dict[int, list[ComponentSite]]] = {}
    for snp_id, delta in phases.items():
        component = components.get(snp_id)
        if component is None:
            continue
        chrom, pos1, ref, alt = parse_snp_id(snp_id)
        site_map.setdefault(chrom, {}).setdefault(pos1 - 1, []).append(
            ComponentSite(ref=ref, alt=alt, delta=int(delta), component=str(component), snp_id=snp_id)
        )
    site_positions = {chrom: sorted(by_pos) for chrom, by_pos in site_map.items()}
    return site_map, site_positions


def load_component_orientations(path: str | Path) -> dict[str, int]:
    out: dict[str, int] = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["component_id"]] = int(row["orientation"])
    return out


def load_component_orientation_table(path: str | Path) -> tuple[dict[str, int], dict[str, str]]:
    orientations: dict[str, int] = {}
    phase_sets: dict[str, str] = {}
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            component_id = row["component_id"]
            orientations[component_id] = int(row["orientation"])
            phase_sets[component_id] = row.get("hic_ps") or component_id
    return orientations, phase_sets


__all__ = [
    "ComponentSite",
    "load_component_orientation_table",
    "load_component_orientations",
    "load_component_site_map",
]
