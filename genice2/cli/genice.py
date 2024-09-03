import sys
import argparse as ap
from collections import defaultdict
from genice2.plugin import safe_import, descriptions

# from genice2 import __version__
from genice2.cli import SmartFormatter, logger_setup, help_water, help_format
from genice2.genice import GenIce
from genice2.valueparser import plugin_option_parser, parse_guest
from genice2.decorators import timeit, banner
import pickle
import numpy as np
from importlib.metadata import version

# 遅延評価。descriptions()関数は重いので、必要なければ呼びたくない。


def help_type():
    return (
        "R|Crystal type (1c, 1h, etc. See https://github.com/vitroid/GenIce for available ice structures.)\n\nIf you want to analyze your own structures, please try analice tool.\n\n"
        + descriptions("lattice", width=55)
    )


def help_guest():
    return (
        "R|Specify guest(s) in the cage type. (D=empty, T=co2*0.5+me*0.3, etc.)\n\n"
        + descriptions("molecule", water=False, width=55)
    )


def getoptions():
    try:
        genice2_version = version("genice2")
    except:
        genice2_version = "2.*.*.*"
    parser = ap.ArgumentParser(
        description=f"GenIce is a swiss army knife to generate hydrogen-disordered ice structures. (version {genice2_version})",
        prog="genice2",
        formatter_class=SmartFormatter,
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version="%(prog)s {0}".format(genice2_version),
    ),
    parser.add_argument(
        "--rep",
        "-r",
        nargs=3,
        type=int,
        dest="rep",
        default=[1, 1, 1],
        help="Repeat the unit cell along a, b, and c axes. [1,1,1]",
    )
    parser.add_argument(
        "--reshape",
        "-R",
        type=str,
        dest="reshape",
        default="",
        help=(
            "Convert the unit cell shape by specifying the new (a,b,c) set from the original (a,b,c) of the unit cell. "
            + "The combination of (a,b,c) is specified by nine integers. "
            + "For example, '--reshape 3,0,0,0,2,0,0,0,1' specifies that the new cell vectors are (3a, 2b, c), which is equivalent to '--rep 3 2 1'."
        ),
    )
    parser.add_argument(
        "--shift",
        "-S",
        nargs=3,
        type=float,
        dest="shift",
        default=[0.0, 0.0, 0.0],
        help="Shift the unit cell along a, b, and c axes. (0.5==half cell) [0,0,0]",
    )
    parser.add_argument(
        "--dens",
        "-d",
        type=float,
        dest="dens",
        default=-1,
        help="Specify the ice density in g/cm3 (Guests are not included.)",
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
        "--seed", "-s", type=int, dest="seed", default=1000, help="Random seed [1000]"
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
        "--water",
        "-w",
        dest="water",
        default="tip3p",
        metavar="model",
        help=help_water(),
    )
    parser.add_argument(
        "--guest",
        "-g",
        # nargs=1,
        dest="guests",
        metavar="14=me",
        action="append",
        help=help_guest(),
    )
    parser.add_argument(
        "--Guest",
        "-G",
        # nargs=1,
        dest="spot_guests",
        metavar="13=me",
        action="append",
        help="Specify guest in the specific cage. (13=me, 32=co2, etc.)",
    )
    parser.add_argument(
        "--Group",
        "-H",
        # nargs=1,
        dest="groups",
        metavar="13=bu-:0",
        action="append",
        help="Specify the group. (-H 13=bu-:0, etc.)",
    )
    parser.add_argument(
        "--anion",
        "-a",
        # nargs=1,
        dest="anions",
        metavar="3=Cl",
        action="append",
        help="Specify a monatomic anion that replaces a water molecule. (3=Cl, 39=F, etc.)",
    )
    parser.add_argument(
        "--cation",
        "-c",
        # nargs=1,
        dest="cations",
        metavar="3=Na",
        action="append",
        help="Specify a monatomic cation that replaces a water molecule. (3=Na, 39=NH4, etc.)",
    )
    # parser.add_argument('--visual',
    #                     dest='visual',
    #                     default="",
    #                     metavar="visual",
    # help='Specify the yaplot file to store the depolarization paths. [""]')
    parser.add_argument(
        "--depol",
        dest="depol",
        default="strict",
        help='Depolarization. (strict, optimal, or none) ["strict"]',
    )
    parser.add_argument(
        "--asis",
        action="store_true",
        dest="asis",
        default=False,
        help="Assumes all given HB pairs to be fixed. No shuffle and no depolarization.",
    )
    parser.add_argument(
        "--debug",
        "-D",
        action="store_true",
        dest="debug",
        help="Output debugging info.",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        dest="quiet",
        help="Do not output progress messages.",
    )
    parser.add_argument(
        "--assess_cages",
        "-A",
        action="store_true",
        dest="assess_cages",
        help="Assess the locations of cages based on the HB network topology. Note: it may fail when the unit cell is too small.",
    )
    parser.add_argument("Type", help=help_type())
    return parser.parse_args()


@timeit
@banner
def main():
    """
    GenIce
    """
    # Module-loading paths
    # 1. Look for the modules in the current working directory
    sys.path.append(".")

    # Parse options
    options = getoptions()

    # Set verbosity level
    logger = logger_setup(options.debug, options.quiet)
    logger.debug("Debug mode.")

    lattice_type = options.Type

    seed = options.seed
    np.random.seed(seed)
    # self.seed = seed  # used in tilecycles
    if options.reshape != "":
        reshape = np.array([int(x) for x in options.reshape.split(",")]).reshape(3, 3)
    else:
        reshape = np.diag(options.rep)
    sh = options.shift
    density = options.dens
    asis = options.asis
    anions = dict()
    if options.anions is not None:
        logger.info(options.anions)
        for v in options.anions:
            key, value = v.split("=")
            anions[int(key)] = value
    cations = dict()
    if options.cations is not None:
        for v in options.cations:
            key, value = v.split("=")
            cations[int(key)] = value
    spot_guests = dict()
    if options.spot_guests is not None:
        for v in options.spot_guests:
            key, value = v.split("=")
            spot_guests[int(key)] = value
    groups = dict()
    if options.groups is not None:
        for v in options.groups:
            key, value = v.split("=")
            groups[int(key)] = value

    lattice_type, lattice_options = plugin_option_parser(options.Type)
    logger.debug("Lattice: {0}".format(lattice_type))
    assert lattice_type is not None

    signature = "Command line: {0}".format(" ".join(sys.argv))

    # Initialize the Lattice class with arguments which are required for
    # plugins.
    lat = GenIce(
        safe_import("lattice", lattice_type).Lattice(**lattice_options),
        signature=signature,
        density=density,
        reshape=reshape,
        cations=cations,
        anions=anions,
        spot_guests=spot_guests,
        spot_groups=groups,
        asis=asis,
        shift=sh,
        rep=None,
    )

    guests = defaultdict(dict)
    if options.guests is not None:
        logger.info(options.guests)
        for guest_spec in options.guests:
            parse_guest(guests, guest_spec)
    noise = options.noise
    depol = options.depol
    assess_cages = options.assess_cages
    file_format, format_options = plugin_option_parser(options.format)

    water_type, water_options = plugin_option_parser(options.water)
    logger.debug("Water type: {0}".format(water_type))
    water = safe_import("molecule", water_type).Molecule(**water_options)

    # Main part of the program is contained in th Formatter object. (See
    # formats/)
    logger.debug("Output file format: {0}".format(file_format))
    formatter_module = safe_import("format", file_format)
    formatter = formatter_module.Format(**format_options)

    del options  # Dispose for safety.

    result = lat.generate_ice(
        water=water,
        guests=guests,
        formatter=formatter,
        noise=noise,
        depol=depol,
        assess_cages=assess_cages,
    )
    if isinstance(result, bytes):
        sys.stdout.buffer.write(result)
    elif isinstance(result, str):
        sys.stdout.write(result)
    elif result is not None:
        pickle.dump(result, sys.stdout.buffer)


if __name__ == "__main__":
    main()
