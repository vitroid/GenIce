import sys

# import logging
import argparse as ap
from genice2.cli import SmartFormatter, help_format, logger_setup
from genice2 import load, analice
from importlib.metadata import version

# from genice2.plugin import descriptions
from genice2.plugin import safe_import, descriptions
from genice2.valueparser import plugin_option_parser
import numpy as np


def help_file():
    return (
        "R|Input file(s). File type is estimated from the suffix. Files of different types cannot be read at a time. File type can be specified explicitly with -s option.\n\n"
        + descriptions("loader")
    )


def getoptions():
    try:
        genice2_version = version("genice2")
    except:
        genice2_version = "2.*.*.*"
    parser = ap.ArgumentParser(
        description=f"GenIce is a swiss army knife to generate hydrogen-disordered ice structures. (version {genice2_version})",
        prog="analice2",
        usage="%(prog)s [options]",
        formatter_class=SmartFormatter,
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version="%(prog)s {0}".format(genice2_version),
    )
    parser.add_argument(
        "--format",
        "-f",
        dest="format",
        default="gromacs",
        metavar="name",
        help=help_format(),
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output",
        metavar="%04d.gro",
        help="Output in separate files.",
    )
    parser.add_argument(
        "--water",
        "-w",
        dest="water",
        default="tip3p",
        metavar="model",
        help="Replace water model. (tip3p, tip4p, etc.) [tip3p]",
    )
    parser.add_argument(
        "--oxygen",
        "-O",
        dest="oatom",
        metavar="OW",
        default="O",
        help='Specify atom name of oxygen in input Gromacs file. ("O")',
    )
    parser.add_argument(
        "--hydrogen",
        "-H",
        dest="hatom",
        metavar="HW[12]",
        default="H",
        help='Specify atom name (regexp) of hydrogen in input Gromacs file. ("H")',
    )
    parser.add_argument(
        "--suffix",
        "-s",
        dest="suffix",
        metavar="gro",
        default=None,
        help="Specify the file suffix explicitly. ((None)",
    )
    parser.add_argument(
        "--filerange",
        dest="filerange",
        metavar="[from:]below[:interval]",
        default="0:1000000",
        help='Specify the number range for the input filename. ("0:1000000")',
    )
    parser.add_argument(
        "--framerange",
        dest="framerange",
        metavar="[from:]below[:interval]",
        default="0:1000000",
        help='Specify the number range for the input frames. ("0:1000000")',
    )
    parser.add_argument(
        "--debug", "-D", action="count", dest="debug", help="Output debugging info."
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        dest="quiet",
        help="Do not output progress messages.",
    )
    parser.add_argument(
        "--seed", type=int, dest="seed", default=1000, help="Random seed [1000]"
    )
    parser.add_argument(
        "--add_noise",
        type=float,
        dest="noise",
        default=0.0,
        metavar="percent",
        help="Add a Gauss noise with given width (SD) to the molecular positions of water. The value 1 corresponds to 1 percent of the molecular diameter of water.",
    )
    parser.add_argument(
        "--avgspan",
        "-v",
        type=float,
        dest="avgspan",
        default=0,
        metavar="1",
        help="Output mean atomic positions of a given time span so as to remove fast librational motions and to make a smooth video. The values 0 and 1 specify no averaging.",
    )
    parser.add_argument("File", help=help_file())
    return parser.parse_args()


def main():
    # Module-loading paths
    # 1. Look for the modules in the current working directory
    sys.path.append(".")

    # Parse options
    options = getoptions()

    logger = logger_setup(options.debug, options.quiet)
    logger.debug("Debug mode.")

    logger.debug(options.File)

    file_format, format_options = plugin_option_parser(options.format)
    logger.debug("Output file format: {0}".format(file_format))
    formatter_module = safe_import("format", file_format)
    formatter = formatter_module.Format(**format_options)

    water_type, water_options = plugin_option_parser(options.water)
    logger.debug("Water type: {0}".format(water_type))
    water = safe_import("molecule", water_type).Molecule(**water_options)

    seed = options.seed
    np.random.seed(seed)

    oname = options.oatom
    hname = options.hatom
    filename = options.File
    noise = options.noise
    avgspan = options.avgspan
    filerange = options.filerange
    framerange = options.framerange
    suffix = options.suffix
    if options.output is None:
        output = None
        stdout = None
    else:
        output = options.output
        stdout = sys.stdout
    signature = "Command line: {0}".format(" ".join(sys.argv))

    logger.debug(filerange)
    logger.debug(framerange)
    logger.debug(oname)
    logger.debug(hname)
    logger.debug(suffix)
    logger.info("Output:{0}".format(output))

    del options  # Dispose for safety.

    for i, (oatoms, hatoms, cellmat) in enumerate(
        load.average(
            lambda: load.iterate(
                filename, oname, hname, filerange, framerange, suffix=suffix
            ),
            span=avgspan,
        )
    ):
        lattice_info = load.make_lattice_info(oatoms, hatoms, cellmat)

        result = analice.AnalIce(lattice_info, signature=signature).analyze_ice(
            water=water,
            formatter=formatter,
            noise=noise,
        )

        if result is not None:
            if output is not None:
                if isinstance(result, bytes):
                    # binary mode
                    # redirect
                    sys.stdout = open(output % i, "wb")  # .buffer
                    logger.debug(f"Binary output to a file: {output%i}")
                else:
                    sys.stdout = open(output % i, "w")
                    logger.debug(f"Text output to a file: {output%i}")
            sys.stdout.write(result)
    if stdout is not None:
        # recover stdout
        sys.stdout = stdout


if __name__ == "__main__":
    main()
