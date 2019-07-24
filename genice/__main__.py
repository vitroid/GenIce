#!/usr/bin/env python3
# -*- python -*-

import os
import sys
import argparse as ap
import logging
from genice.importer import safe_import
from genice import lattice, __version__, load
import random
import numpy as np


def getoptions():
    parser = ap.ArgumentParser(description='GenIce is a swiss army knife to generate hydrogen-disordered ice structures. (version {0})'.format(__version__), prog='genice')
    parser.add_argument('--version',
                        '-V',
                        action='version',
                        version='%(prog)s {0}'.format(__version__))
    parser.add_argument('--rep',
                        '-r',
                        nargs=3,
                        type=int,
                        dest='rep',
                        default=[1, 1, 1],
                        help='Repeat the unit cell in x,y, and z directions. [1,1,1]')
    parser.add_argument('--dens',
                        '-d',
                        type=float,
                        dest='dens',
                        default=-1,
                        help='Specify the ice density in g/cm3')
    parser.add_argument('--add_noise',
                        type=float,
                        dest='noise',
                        default=0.,
                        metavar='percent',
                        help='Add a Gauss noise with given width (SD) to the molecular positions of water. The value 1 corresponds to 1 percent of the molecular diameter of water.')
    parser.add_argument('--seed',
                        '-s',
                        type=int,
                        dest='seed',
                        default=1000,
                        help='Random seed [1000]')
    parser.add_argument('--format',
                        '-f',
                        dest='format',
                        default="gromacs",
                        metavar="gmeqdypoc",
                        help='Specify file format [g(romacs)|m(dview)|e(uler)|q(uaternion)|d(igraph)|y(aplot)|p(ython module)|o(penScad)|c(entersofmass)|r(elative com)] [gromacs]')
    parser.add_argument('--water',
                        '-w',
                        dest='water',
                        default="tip3p",
                        metavar="model",
                        help='Specify water model. (tip3p, tip4p, etc.) [tip3p]')
    parser.add_argument('--guest',
                        '-g',
                        nargs=1,
                        dest='guests',
                        metavar="D=empty",
                        action="append",
                        help='Specify guest(s) in the cage type. (D=empty, T=co2*0.5+me*0.3, etc.)')
    parser.add_argument('--Guest',
                        '-G',
                        nargs=1,
                        dest='spot_guests',
                        metavar="13=me",
                        action="append",
                        help='Specify guest in the specific cage. (13=me, 32=co2, etc.)')
    parser.add_argument('--Group',
                        '-H',
                        nargs=1,
                        dest='groups',
                        metavar="13=bu-:0",
                        action="append",
                        help='Specify the group. (-H 13=bu-:0, etc.)')
    parser.add_argument('--anion',
                        '-a',
                        nargs=1,
                        dest='anions',
                        metavar="3=Cl",
                        action="append",
                        help='Specify a monatomic anion that replaces a water molecule. (3=Cl, 39=F, etc.)')
    parser.add_argument('--cation',
                        '-c',
                        nargs=1,
                        dest='cations',
                        metavar="3=Na",
                        action="append",
                        help='Specify a monatomic cation that replaces a water molecule. (3=Na, 39=NH4, etc.)')
    parser.add_argument('--visual',
                        dest='visual',
                        default="",
                        metavar="visual",
                        help='Specify the yaplot file to store the depolarization paths. [""]')
    parser.add_argument('--nodep',
                        action='store_true',
                        dest='nodep',
                        default=False,
                        help='No depolarization.')
    parser.add_argument('--asis',
                        action='store_true',
                        dest='asis',
                        default=False,
                        help='Assumes all given HB pairs to be fixed. No shuffle and no depolarization.')
    parser.add_argument('--debug',
                        '-D',
                        action='store_true',
                        dest='debug',
                        help='Output debugging info.')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        dest='quiet',
                        help='Do not output progress messages.')
    parser.add_argument('Type',
                        help='Crystal type (1c,1h,etc. See https://github.com/vitroid/GenIce for available ice structures.)')
    return parser.parse_args()


def getoptions_analice():
    parser = ap.ArgumentParser(description='GenIce is a swiss army knife to generate hydrogen-disordered ice structures. (version {0})'.format(__version__), prog='analice', usage='%(prog)s [options]')
    parser.add_argument('--version',
                        '-V',
                        action='version',
                        version='%(prog)s {0}'.format(__version__))
    parser.add_argument('--format',
                        '-f',
                        dest='format',
                        default="gromacs",
                        metavar="gmeqdypoc",
                        help='Specify file format [g(romacs)|m(dview)|e(uler)|q(uaternion)|d(igraph)|y(aplot)|p(ython module)|o(penScad)|c(entersofmass)|r(elative com)] [gromacs]')
    parser.add_argument('--output',
                        '-o',
                        dest='output',
                        metavar="%04d.gro",
                        help='Output in separate files.')
    parser.add_argument('--water',
                        '-w',
                        dest='water',
                        default="tip3p",
                        metavar="model",
                        help='Replace water model. (tip3p, tip4p, etc.) [tip3p]')
    parser.add_argument('--oxygen',
                        '-O',
                        dest='oatom',
                        metavar="OW",
                        default="O",
                        help='Specify atom name of oxygen in input Gromacs file. ("O")')
    parser.add_argument('--hydrogen',
                        '-H',
                        dest='hatom',
                        metavar="HW[12]",
                        default="H",
                        help='Specify atom name (regexp) of hydrogen in input Gromacs file. ("H")')
    parser.add_argument('--suffix',
                        '-s',
                        dest='suffix',
                        metavar="gro",
                        default=None,
                        help='Override the file suffix. (None)')
    parser.add_argument('--filerange',
                        dest='filerange',
                        metavar="[from:]below[:interval]",
                        default="0:1000000",
                        help='Specify the number range for the input filename. ("0:1000000")')
    parser.add_argument('--framerange',
                        dest='framerange',
                        metavar="[from:]below[:interval]",
                        default="0:1000000",
                        help='Specify the number range for the input frames. ("0:1000000")')
    parser.add_argument('--debug',
                        '-D',
                        action='count',
                        dest='debug',
                        help='Output debugging info.')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        dest='quiet',
                        help='Do not output progress messages.')
    parser.add_argument('--add_noise',
                        type=float,
                        dest='noise',
                        default=0.,
                        metavar='percent',
                        help='Add a Gauss noise with given width (SD) to the molecular positions of water. The value 1 corresponds to 1 percent of the molecular diameter of water.')
    parser.add_argument('--avgspan', '-v',
                        type=float,
                        dest='avgspan',
                        default=0,
                        metavar='1',
                        help='Average atomic positions in water molecules so as to remove fast librational motions and to make a smooth video. Specify the average span. The values 0 and 1 specify no averaging.')
    parser.add_argument('File',
                        help='Gromacs file.')
    return parser.parse_args()


def main():
    # Module-loading paths
    # 1. Look for the modules in the current working directory
    sys.path.append(".")

    # Parse options
    if sys.argv[0].find("analice") >= 0:
        options = getoptions_analice()
        mode = "analice"
    else:
        options = getoptions()
        mode = "genice"

    # Set verbosity level
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    elif options.quiet:
        logging.basicConfig(level=logging.WARN,
                            format="%(levelname)s %(message)s")
    else:
        # normal
        logging.basicConfig(level=logging.INFO,
                            format="%(levelname)s %(message)s")
    logger = logging.getLogger()
    logger.debug("Debug mode.")

    if mode == "genice":
        logger.debug(options.Type)

        lattice_type = options.Type
        seed = options.seed
        rep = options.rep
        density = options.dens
        asis = options.asis
        anions = dict()
        if options.anions is not None:
            logger.info(options.anions)
            for v in options.anions:
                key, value = v[0].split("=")
                anions[int(key)] = value
        cations = dict()
        if options.cations is not None:
            for v in options.cations:
                key, value = v[0].split("=")
                cations[int(key)] = value
        spot_guests = dict()
        if options.spot_guests is not None:
            for v in options.spot_guests:
                key, value = v[0].split("=")
                spot_guests[int(key)] = value
        groups = dict()
        if options.groups is not None:
            for v in options.groups:
                key, value = v[0].split("=")
                groups[int(key)] = value

        # Set random seeds
        random.seed(seed)
        np.random.seed(seed)

        logger.debug("Lattice: {0}".format(lattice_type))
        assert lattice_type is not None

        # Initialize the Lattice class with arguments which are required for plugins.
        lat = lattice.Lattice(safe_import("lattice", lattice_type),
                              density=density,
                              rep=rep,
                              cations=cations,
                              anions=anions,
                              spot_guests=spot_guests,
                              spot_groups=groups,
                              asis=asis,
                              )

        water_type = options.water
        guests = options.guests
        noise = options.noise
        depolarize = not options.nodep
        file_format = options.format

        # Main part of the program is contained in th Formatter object. (See formats/)
        logger.debug("Output file format: {0}".format(file_format))
        formatter = safe_import("format", file_format)

        if options.visual != "":
            record_depolarization_path = open(options.visual, "w")
        else:
            record_depolarization_path = None

        del options  # Dispose for safety.

        lat.generate_ice(water_type=water_type,
                         guests=guests,
                         formatter=formatter,
                         record_depolarization_path=record_depolarization_path,
                         noise=noise,
                         depolarize=depolarize,
                         )
    else:  # analice
        logger.debug(options.File)

        water_type = options.water
        file_format = options.format
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

        logger.debug(filerange)
        logger.debug(framerange)
        logger.debug(oname)
        logger.debug(hname)
        logger.debug(suffix)
        logger.info("Output:{0}".format(output))

        del options  # Dispose for safety.

        for i, w in enumerate(load.iterate(filename, oname, hname, filerange, framerange, suffix=suffix, avgspan=avgspan)):
            # Main part of the program is contained in th Formatter object. (See formats/)
            logger.debug("Output file format: {0}".format(file_format))
            formatter = safe_import("format", file_format)
            lat = lattice.Lattice(w)
            if output is not None:
                sys.stdout = open(output % i, "w")
            lat.analyze_ice(water_type=water_type,
                            formatter=formatter,
                            noise=noise,
                            )
        if stdout is not None:
            # recover stdout
            sys.stdout = stdout


if __name__ == "__main__":
    main()
