import pytest

from phirefly.pipeline import build_parser


def test_short_help_keeps_common_options_and_full_help_is_reversible(capsys):
    parser = build_parser()
    before = parser.format_help()
    assert '--bam' in before and '--threads' in before and '--help-all' in before
    assert '--risk-threshold' not in before and '--truth-vcf' not in before
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(['--help-all'])
    assert exc.value.code == 0
    full = capsys.readouterr().out
    assert '--risk-threshold' in full and '--truth-vcf' in full
    assert 'not a calibrated error probability' in full
    assert parser.format_help() == before


def test_hidden_options_still_parse_and_defaults_do_not_change():
    parser = build_parser()
    required = ['--bam', 'reads.bam', '--vcf', 'sites.vcf.gz',
                '--region', 'chr1:1-100', '--out-dir', 'out']
    args = parser.parse_args(required)
    assert (args.algorithm, args.batch_size, args.n_iter, args.seed) == ('CFC', 100, 1000, 23)
    assert (args.risk_threshold, args.soft_risk_threshold, args.soft_edge_scale) == (.45, .75, .25)
    changed = parser.parse_args(required + ['--risk-threshold', '0.4', '--n-iter', '2000'])
    assert changed.risk_threshold == .4 and changed.n_iter == 2000
    assert all(action.help for action in parser._actions)
