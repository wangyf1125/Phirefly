import sys

import pytest

from phirefly.core import runtime
from phirefly.pipeline import build_parser, run_pipeline


def test_foreign_python_path_is_rejected(monkeypatch):
    monkeypatch.setattr(sys, 'path', ['/apps/lib/python3.9/site-packages'])
    with pytest.raises(ValueError, match='Foreign Python packages'):
        runtime.check_python_environment()


def test_missing_tools_fail_before_output_creation(tmp_path, monkeypatch):
    monkeypatch.setattr(runtime.shutil, 'which', lambda _: None)
    args = build_parser().parse_args([
        '--bam', 'unused.bam', '--vcf', 'unused.vcf', '--region', 'chr1:1-100',
        '--out-dir', str(tmp_path / 'out')])
    with pytest.raises(ValueError, match='Missing executable.*bgzip.*tabix'):
        run_pipeline(args)
    assert not (tmp_path / 'out').exists()


def test_runtime_checks_requested_tools_and_dependency(monkeypatch):
    monkeypatch.setattr(runtime.shutil, 'which', lambda command: '/tools/' + command)
    monkeypatch.setattr(runtime, 'find_spec', lambda _: None)
    with pytest.raises(ValueError, match='MindQuantum is not installed'):
        runtime.check_runtime('custom-bgzip', 'custom-tabix')
    monkeypatch.setattr(runtime, 'find_spec', lambda _: object())
    assert runtime.check_runtime('custom-bgzip', 'custom-tabix') == {
        'bgzip': '/tools/custom-bgzip', 'tabix': '/tools/custom-tabix'}
