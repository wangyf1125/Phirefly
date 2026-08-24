#!/usr/bin/env python3
"""Export SNP phase predictions as a phased VCF."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pysam

try:
    from .extract.store import iter_observation_pairs
except ImportError:  # pragma: no cover - direct script execution fallback
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from phirefly.extract.store import iter_observation_pairs  # type: ignore


def load_phases(path: str | Path) -> dict[str, int]:
    phases, _phase_sets = load_phase_records(path)
    return phases


def load_phase_records(path: str | Path) -> tuple[dict[str, int], dict[str, int]]:
    phases: dict[str, int] = {}
    phase_sets: dict[str, int] = {}
    with Path(path).open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {key: i for i, key in enumerate(header)}
        ps_col = idx.get("ps", idx.get("PS"))
        for line in f:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            snp_id = row[idx["snp_id"]]
            phases[snp_id] = int(row[idx["delta"]])
            if ps_col is not None and row[ps_col] not in {"", "NA", "."}:
                phase_sets[snp_id] = int(float(row[ps_col]))
    return phases, phase_sets


def snp_pos_from_id(snp_id: str) -> int:
    return int(snp_id.split(":", 2)[1])


def load_component_phase_sets(observations_path: str | Path) -> dict[str, int]:
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != x:
            next_x = parent[x]
            parent[x] = root
            x = next_x
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for read_id, snp_id in iter_observation_pairs(str(observations_path)):
        union("R:" + read_id, "S:" + snp_id)

    snp_to_component: dict[str, str] = {}
    component_min_pos: dict[str, int] = {}
    for node in list(parent):
        if not node.startswith("S:"):
            continue
        snp_id = node[2:]
        component = find(node)
        snp_to_component[snp_id] = component
        pos = snp_pos_from_id(snp_id)
        component_min_pos[component] = min(component_min_pos.get(component, pos), pos)
    return {snp_id: component_min_pos[component] for snp_id, component in snp_to_component.items()}


def parse_region(region: str) -> tuple[str, int, int]:
    chrom, rest = region.split(":", 1)
    start_s, end_s = rest.replace(",", "").split("-", 1)
    return chrom, int(start_s) - 1, int(end_s)


def region_records(vcf: pysam.VariantFile, chrom: str, start0: int, end0: int):
    try:
        yield from vcf.fetch(chrom, start0, end0)
    except (ValueError, OSError):
        for rec in vcf:
            if rec.chrom == chrom and start0 < rec.pos <= end0:
                yield rec


def export_phased_vcf(
    input_vcf: str | Path,
    pred_snp_phases: str | Path,
    region: str,
    output_vcf: str | Path,
    sample: str | None = None,
    ps_start: int | None = None,
    observations: str | Path | None = None,
) -> tuple[int, int, Path]:
    """Copy an input VCF header and write phased GT/PS calls for predicted SNPs."""

    phases, phase_sets = load_phase_records(pred_snp_phases)
    ps_by_snp = load_component_phase_sets(observations) if observations else {}
    chrom, start0, end0 = parse_region(region)

    invcf = pysam.VariantFile(str(input_vcf))
    samples = list(invcf.header.samples)
    if not samples:
        raise SystemExit(f"input VCF has no samples: {input_vcf}")
    sample_name = sample or samples[0]
    if sample_name not in invcf.header.samples:
        raise SystemExit(f"sample {sample_name} not found in {input_vcf}")

    invcf.subset_samples([sample_name])
    out_header = invcf.header.copy()
    if "PS" not in out_header.formats:
        out_header.formats.add("PS", 1, "Integer", "Phase set")
    out_path = Path(output_vcf)
    out_mode = "wz" if str(out_path).endswith(".gz") else "w"
    default_ps = ps_start if ps_start is not None else start0 + 1

    written = 0
    skipped = 0
    with pysam.VariantFile(str(out_path), out_mode, header=out_header) as out:
        for rec in region_records(invcf, chrom, start0, end0):
            if len(rec.ref) != 1 or len(rec.alts or []) != 1 or len(rec.alts[0]) != 1:
                skipped += 1
                continue
            snp_id = f"{rec.chrom}:{rec.pos}:{rec.ref}>{rec.alts[0]}"
            delta = phases.get(snp_id)
            if delta is None:
                skipped += 1
                continue
            if delta == 1:
                gt = (0, 1)
            elif delta == -1:
                gt = (1, 0)
            else:
                raise SystemExit(f"invalid delta for {snp_id}: {delta}")
            out_rec = rec.copy()
            out_rec.translate(out_header)
            out_rec.samples[sample_name]["GT"] = gt
            out_rec.samples[sample_name].phased = True
            out_rec.samples[sample_name]["PS"] = int(phase_sets.get(snp_id, ps_by_snp.get(snp_id, default_ps)))
            out.write(out_rec)
            written += 1
    invcf.close()
    return written, skipped, out_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-vcf", required=True)
    parser.add_argument("--sample", default=None, help="Sample to write. Default: first sample in input VCF.")
    parser.add_argument("--pred-snp-phases", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--ps-start", type=int, default=None)
    parser.add_argument("--observations", default=None)
    args = parser.parse_args()

    written, skipped, out_path = export_phased_vcf(
        input_vcf=args.input_vcf,
        pred_snp_phases=args.pred_snp_phases,
        region=args.region,
        output_vcf=args.output_vcf,
        sample=args.sample,
        ps_start=args.ps_start,
        observations=args.observations,
    )
    print(f"written\t{written}")
    print(f"skipped\t{skipped}")
    print(f"output_vcf\t{out_path}")


if __name__ == "__main__":
    main()


__all__ = [
    "export_phased_vcf",
    "load_component_phase_sets",
    "load_phase_records",
    "load_phases",
    "parse_region",
    "region_records",
]
