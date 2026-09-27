import sys
from types import SimpleNamespace

import pytest

from phirefly import __version__, cli


def test_version_does_not_import_commands(monkeypatch, capsys):
    def unexpected(*args):
        raise AssertionError("version should not import command modules")
    monkeypatch.setattr(cli, "import_module", unexpected)
    monkeypatch.setattr(sys, "argv", ["phirefly", "--version"])
    with pytest.raises(SystemExit) as exc:
        cli.main()
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip() == f"phirefly {__version__}"


def test_dispatch_restores_argv_even_on_error(monkeypatch):
    argv = ["phirefly", "run", "--help"]
    def command():
        assert sys.argv == ["phirefly run", "--help"]
        raise SystemExit(0)
    monkeypatch.setattr(cli, "import_module", lambda name: SimpleNamespace(main=command))
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(SystemExit):
        cli.main()
    assert sys.argv is argv
