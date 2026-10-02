"""Short default help without removing expert command-line options."""

import argparse


class FullHelp(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        hidden = parser._advanced_help
        try:
            for action, help_text in hidden:
                action.help = help_text
            parser.print_help()
        finally:
            for action, _ in hidden:
                action.help = argparse.SUPPRESS
        parser.exit()


def compact_help(parser, common):
    parser.formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser.add_argument('--help-all', action=FullHelp, nargs=0,
                        help='Show all options, including advanced settings, and exit.')
    parser._advanced_help = []
    for action in parser._actions:
        if action.dest not in {*common, 'help', 'help_all'}:
            parser._advanced_help.append((action, action.help))
            action.help = argparse.SUPPRESS
    return parser
