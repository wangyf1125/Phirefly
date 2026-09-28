"""Umbrella command for Phirefly public CLI entry points."""

from __future__ import annotations

import argparse
import sys
from importlib import import_module

from . import __version__
from .core.runtime import check_python_environment

COMMANDS = {
    "run": "pipeline",
    "extract": "extract.cli",
    "msf": "msf",
    "phaselet-qaia": "phaselet.cli",
    "vcf": "vcf",
    "benchmark": "metrics",
    "hic-edges": "hic.component_edges",
    "hic-orient": "hic.orient",
    "hic-apply": "hic.apply",
    "hic-qc": "hic.qc",
}


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="phirefly",
        description="Phirefly phaselet-hyperread phasing command group.",
    )
    parser.add_argument("--version", action="version", version=f"phirefly {__version__}")
    parser.add_argument("command", nargs="?", choices=sorted(COMMANDS))
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        return
    previous = sys.argv
    try:
        check_python_environment()
        command = import_module(f"phirefly.{COMMANDS[args.command]}")
        sys.argv = [f"phirefly {args.command}", *args.args]
        command.main()
    except ValueError as error:
        parser.exit(2, f"phirefly: {error}\n")
    finally:
        sys.argv = previous


if __name__ == "__main__":
    main()
