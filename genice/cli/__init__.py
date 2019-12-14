

import argparse as ap
from textwrap import wrap
from genice.plugin import descriptions

class SmartFormatter(ap.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            # return text[2:].splitlines()
            return [line for L in text[2:].splitlines() for line in wrap(L, width=55, drop_whitespace=False) ]
        # this is the RawTextHelpFormatter._split_lines
        return ap.HelpFormatter._split_lines(self, text, width)
    def _get_help_string(self, action):
        if callable(action.help):
            return action.help()
        return action.help

def help_format():
    return 'R|Specify the output file format. [gromacs]\n\n'+descriptions("format", width=55)
