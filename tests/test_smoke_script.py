import os
from pathlib import Path
import shutil
import subprocess

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / 'examples/mhc_smoke/run_mhc_smoke.sh'


@pytest.mark.parametrize('selected', ['path', 'conda', 'venv', 'explicit', 'missing'])
def test_smoke_interpreter_selection(tmp_path, selected):
    bash = shutil.which('bash')
    if bash is None:
        pytest.skip('Smoke script requires bash')
    env = os.environ.copy()
    for name in ('PHIREFLY_PYTHON', 'CONDA_PREFIX', 'VIRTUAL_ENV', 'BASH_ENV', 'ENV'):
        env.pop(name, None)
    marker = tmp_path / 'selected-python'
    env['SELECTED_PYTHON'] = str(marker)
    interpreters = {}
    for name in ('path', 'conda', 'venv', 'explicit'):
        executable = tmp_path / (name + ' env') / 'bin/python'
        executable.parent.mkdir(parents=True)
        executable.write_text(
            '#!/bin/sh\nprintf "%s\\n" "$0" > "$SELECTED_PYTHON"\nexit 42\n')
        executable.chmod(0o755)
        interpreters[name] = executable
    # Simulate a startup file putting an unrelated Python before the active one.
    startup = tmp_path / 'bash-env'
    startup.write_text('export PATH="$WRONG_PYTHON_BIN:$PATH"\n')
    env['BASH_ENV'] = str(startup)
    env['WRONG_PYTHON_BIN'] = str(interpreters['path'].parent)
    if selected != 'path':
        env['CONDA_PREFIX'] = str(interpreters['conda'].parent.parent)
    if selected in ('venv', 'explicit', 'missing'):
        env['VIRTUAL_ENV'] = str(interpreters['venv'].parent.parent)
    if selected in ('explicit', 'missing'):
        env['PHIREFLY_PYTHON'] = str(
            interpreters['explicit'] if selected == 'explicit' else tmp_path / 'missing')
    result = subprocess.run([bash, str(SCRIPT)], env=env, text=True, capture_output=True)
    assert result.returncode == 1
    assert 'Cannot import Phirefly using:' in result.stderr
    if selected == 'missing':
        assert not marker.exists()
        assert env['PHIREFLY_PYTHON'] in result.stderr
    else:
        assert marker.read_text().strip() == str(interpreters[selected])
        if selected != 'path':
            assert str(interpreters[selected]) in result.stderr
