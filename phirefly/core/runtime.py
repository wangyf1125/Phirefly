"""Fail before creating outputs when the runtime environment is unusable."""

from importlib.util import find_spec
import re
import shutil
import sys


def check_python_environment() -> None:
    if sys.version_info[:2] not in {(3, 10), (3, 11)}:
        raise ValueError("Phirefly requires Python 3.10 or 3.11")
    for path in sys.path:
        match = re.search(r"/python(\d+)\.(\d+)/(?:site|dist)-packages", str(path))
        if match and tuple(map(int, match.groups())) != sys.version_info[:2]:
            raise ValueError(
                f"Foreign Python packages on sys.path: {path}. "
                "Unset PYTHONPATH/PYTHONHOME and activate a clean environment; "
                "do not install with sudo or --user.")


def check_runtime(bgzip: str, tabix: str) -> dict[str, str]:
    check_python_environment()
    tools = {name: shutil.which(command) for name, command in
             (("bgzip", bgzip), ("tabix", tabix))}
    missing = [name for name, executable in tools.items() if executable is None]
    if missing:
        raise ValueError(
            "Missing executable(s): " + ", ".join(missing) + ". "
            "Install htslib in the active environment and check PATH, "
            "or specify --bgzip/--tabix with executable paths. No phasing was started.")
    if find_spec("mindquantum") is None:
        raise ValueError(
            "MindQuantum is not installed in this Python environment. "
            "Install Phirefly with the qaia extra or use the full Conda package.")
    return tools
