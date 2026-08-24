"""Umbrella command for Phirefly public CLI entry points."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable


COMMANDS: dict[str, tuple[str, Callable[[], None]]] = {}


def register(name: str, description: str, main: Callable[[], None]) -> None:
    COMMANDS[name] = (description, main)


def load_commands() -> None:
    if COMMANDS:
        return
    from . import metrics, msf, pipeline, vcf
    from .extract import cli as extract_cli
    from .phaselet import cli as phaselet_cli
    from .hic import apply as hic_apply
    from .hic import component_edges as hic_edges
    from .hic import orient as hic_orient
    from .hic import qc as hic_qc

    register("extract", "BAM/VCF to read-SNP observations", extract_cli.main)
    register("msf", "build MSF SNP backbone", msf.main)
    register("run", "BAM/VCF to phased VCF pipeline", pipeline.main)
    register("phaselet-qaia", "phaselet-hyperread QAIA phasing", phaselet_cli.main)
    register("vcf", "export phased VCF from SNP phases", vcf.main)
    register("benchmark", "truth-based benchmark metrics", metrics.main)
    register("hic-edges", "build Hi-C/Pore-C component edges", hic_edges.main)
    register("hic-orient", "orient components from Hi-C/Pore-C edges using MSF", hic_orient.main)
    register("hic-apply", "apply Hi-C/Pore-C component orientations to VCF", hic_apply.main)
    register("hic-qc", "truth-independent Hi-C/Pore-C edge QC", hic_qc.main)


def main() -> None:
    load_commands()
    parser = argparse.ArgumentParser(
        prog="phirefly",
        description="Phirefly phaselet-hyperread phasing command group.",
    )
    parser.add_argument("command", nargs="?", choices=sorted(COMMANDS))
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        return
    _description, command_main = COMMANDS[args.command]
    sys.argv = [f"phirefly {args.command}", *args.args]
    command_main()


if __name__ == "__main__":
    main()
