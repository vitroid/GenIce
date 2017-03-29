#!/usr/bin/env python3
# -*- python -*-

import os
import sys
import argparse  as ap
import logging
from genice.importer import safe_import


def getoptions():
    parser = ap.ArgumentParser(description='')
    parser.add_argument('--rep',  '-r', nargs = 3, type=int,   dest='rep',  default=[2,2,2],
                        help='Repeat the unit cell in x,y, and z directions. [2,2,2]')
    parser.add_argument('--dens', '-d', nargs = 1, type=float, dest='dens', default=(-1,),
                        help='Specify the ice density in g/cm3')
    parser.add_argument('--seed', '-s', nargs = 1, type=int,   dest='seed', default=(1000,),
                        help='Random seed [1000]')
    parser.add_argument('--format', '-f', nargs = 1,           dest='format', default=("gromacs",), metavar="gmeqdypoc",
                        help='Specify file format [g(romacs)|m(dview)|e(uler)|q(uaternion)|d(igraph)|y(aplot)|p(ython module)|o(penScad)|c(entersofmass)|r(elative com)]')
    parser.add_argument('--water', '-w', nargs = 1,           dest='water', default=("tip3p",), metavar="model",
                        help='Specify water model. (tip3p, tip4p, etc.)')
    parser.add_argument('--guest', '-g', nargs = 1,           dest='guests', metavar="D=empty", action="append", 
                        help='Specify guest in the cage. (D=empty, T=co2, etc.)')
    parser.add_argument('--nodep', action='store_true', dest='nodep',
                        help='No depolarization.')
    parser.add_argument('--debug', '-D', action='store_true', dest='debug',
                        help='Output debugging info.')
    parser.add_argument('--quiet', '-q', action='store_true', dest='quiet',
                        help='Do not output progress messages.')
    parser.add_argument('Type', nargs=1,
                       help='Crystal type (1c,1h,etc. See https://github.com/vitroid/GenIce for available ice structures.)')
    return parser.parse_args()


        

def main():
    #prepare user's workarea
    home = os.path.expanduser("~")
    if os.path.exists(home+"/Library/Application Support"): #MacOS
        homegenice = home+"/Library/Application Support/GenIce"
    else:
        homegenice = os.path.expanduser(home + "/.genice") #Other unix
    sys.path.append(homegenice)
    try:
        os.makedirs(homegenice+"/formats")
        os.makedirs(homegenice+"/lattices")
        os.makedirs(homegenice+"/molecules")
    except:
        pass #just ignore when failed.


    #Parse options
    options = getoptions()

    
    #Set verbosity level
    if options.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format="%(asctime)s %(levelname)s %(message)s")
    elif options.quiet:
        logging.basicConfig(level=logging.WARN,
                            format="%(asctime)s %(levelname)s %(message)s")
    else:
        #normal
        logging.basicConfig(level=logging.INFO,
                            format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger()
    logger.debug("Debug mode.")
    logger.debug(options.Type)

    
    #Main part of the program is contained in th Formatter object. (See formats/)
    logger.debug("Format: {0}".format(options.format[0]))
    formatter = safe_import("format", options.format[0])
    f = formatter.Formatter(options)
    f.run(options)
    
    
if __name__ == "__main__":
    main()
