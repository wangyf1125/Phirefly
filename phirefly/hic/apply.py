#!/usr/bin/env python3
"""Apply Hi-C/Pore-C component orientations to a Phirefly phased VCF."""

from __future__ import annotations

import argparse
from pathlib import Path

from ..vcf import export_phased_vcf, load_component_phase_sets, load_phases
from .components import load_component_orientation_table


def write_oriented_snp_phases(
    path: Path,
    snp_phases: str | Path,
    observations: str | Path,
    component_orientations: str | Path,
) -> int:
    phases = load_phases(snp_phases)
    components = load_component_phase_sets(observations)
    orientations, hic_phase_sets = load_component_orientation_table(component_orientations)
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w") as out:
        out.write("snp_id\tdelta\tps\n")
        for snp_id in sorted(phases, key=lambda s: (s.split(":", 1)[0], int(s.split(":", 2)[1]), s)):
            component = str(components.get(snp_id, ""))
            orient = int(orientations.get(component, 1))
            ps = hic_phase_sets.get(component, component)
            out.write(f"{snp_id}\t{int(phases[snp_id]) * orient}\t{ps}\n")
            written += 1
    return written


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--snp-phases", required=True)
    parser.add_argument("--observations", required=True)
    parser.add_argument("--component-orientations", required=True)
    parser.add_argument("--component-edges", default=None, help="Accepted for manifest compatibility; not needed for applying orientations.")
    parser.add_argument("--min-contacts", type=int, default=20)
    parser.add_argument("--min-confidence", type=float, default=0.8)
    parser.add_argument("--min-abs-score", type=float, default=0.0)
    parser.add_argument("--region", required=True)
    parser.add_argument("--output-vcf", required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    output_vcf = Path(args.output_vcf)
    adjusted_phases = output_vcf.with_suffix(".snp_phases.tsv")
    n_snps = write_oriented_snp_phases(
        adjusted_phases,
        snp_phases=args.snp_phases,
        observations=args.observations,
        component_orientations=args.component_orientations,
    )
    written, skipped, out_path = export_phased_vcf(
        input_vcf=args.input_vcf,
        pred_snp_phases=adjusted_phases,
        region=args.region,
        output_vcf=output_vcf,
        sample=args.sample,
        ps_start=1,
        observations=args.observations,
    )
    print(f"oriented_snps\t{n_snps}")
    print(f"written\t{written}")
    print(f"skipped\t{skipped}")
    print(f"output_vcf\t{out_path}")


if __name__ == "__main__":
    main()


__all__ = ["main", "write_oriented_snp_phases"]
