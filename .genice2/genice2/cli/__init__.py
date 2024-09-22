import argparse as ap
from textwrap import wrap
from genice2.plugin import descriptions
import logging
from genice2.decorators import timeit, banner

# common helps for genice and analice


def help_format():
    return "R|Specify the output file format. [gromacs]\n\n" + descriptions(
        "format", width=55
    )


def help_water():
    return "R|Specify the water model. [tip3p]\n\n" + descriptions(
        "molecule", water=True, width=55
    )


class SmartFormatter(ap.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith("R|"):
            # return text[2:].splitlines()
            return [
                line
                for L in text[2:].splitlines()
                for line in wrap(L, width=55, drop_whitespace=False)
            ]
        # this is the RawTextHelpFormatter._split_lines
        return ap.HelpFormatter._split_lines(self, text, width)

    @timeit
    @banner
    def _get_help_string(self, action):
        "General help handler"
        if callable(action.help):
            return action.help()
        return action.help


def logger_setup(debug=False, quiet=False):
    # Set verbosity level
    if debug:
        logging.basicConfig(
            level=logging.DEBUG, format="%(asctime)s %(levelname)s:%(message)s"
        )
    elif quiet:
        logging.basicConfig(level=logging.WARN, format="%(levelname)s:%(message)s")
    else:
        # normal
        logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
    return logging.getLogger()
