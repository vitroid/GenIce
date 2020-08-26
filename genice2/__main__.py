#!/usr/bin/env python3
# -*- python -*-

import sys
import logging
from genice2.importer import safe_import
from genice2 import genice, analice, __version__, load
from genice2.tool import plugin_option_parser
import numpy as np






def main():
    # Module-loading paths
    # 1. Look for the modules in the current working directory
    sys.path.append(".")

    # Parse options
    if sys.argv[0].find("analice") >= 0:
        options = analice.getoptions()
        mode = "analice"
    else:
        options = genice.getoptions()
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

        lattice_type, lattice_options = plugin_option_parser(options.Type)
        logger.debug("Lattice: {0}".format(lattice_type))
        assert lattice_type is not None

        comment = "Command line: {0}".format(" ".join(sys.argv))

        # Initialize the Lattice class with arguments which are required for plugins.
        lat = genice.GenIce(safe_import("lattice", lattice_type).Lattice(**lattice_options),
                            density=density,
                            rep=rep,
                            cations=cations,
                            anions=anions,
                            spot_guests=spot_guests,
                            spot_groups=groups,
                            asis=asis,
                            comment=comment,
        )

        water_type = options.water
        guests = options.guests
        noise = options.noise
        depolarize = not options.nodep
        file_format, format_options = plugin_option_parser(options.format)

        logger.debug("Water type: {0}".format(water_type))
        water = safe_import("molecule", water_type).Molecule()
        # Main part of the program is contained in th Formatter object. (See formats/)
        logger.debug("Output file format: {0}".format(file_format))
        formatter_module = safe_import("format", file_format)
        formatter = formatter_module.Format(**format_options)

        if options.visual != "":
            record_depolarization_path = open(options.visual, "w")
        else:
            record_depolarization_path = None

        del options  # Dispose for safety.

        ice = lat.generate_ice(water=water,
                         guests=guests,
                         formatter=formatter,
                         record_depolarization_path=record_depolarization_path,
                         noise=noise,
                         depolarize=depolarize,
                         seed=seed
                         )
        # この方法の場合、データをそのままうけとれるので、jupyterなんかには適しているが、
        # vpythonの場合にはどうしたらいいんだろう。
        # 何も返却しない、ice==Noneというのが正しい動作か。
        if type(ice) is bytes:
            sys.stdout.buffer.write(ice)
        else:
            sys.stdout.write(ice)

    else:  # analice
        logger.debug(options.File)

        water_type = options.water
        file_format, format_options = plugin_option_parser(options.format)
        logger.debug("Output file format: {0}".format(file_format))
        formatter_module = safe_import("format", file_format)
        formatter = formatter_module.Format(**format_options)
        logger.debug("Water type: {0}".format(water_type))
        water_module = safe_import("molecule", water_type)

        oname = options.oatom
        hname = options.hatom
        filename = options.File
        noise = options.noise
        avgspan = options.avgspan
        filerange = options.filerange
        framerange = options.framerange
        suffix = options.suffix
        comment = "Command line: {0}".format(" ".join(sys.argv))
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

        for i, (oatoms, hatoms, cellmat) in enumerate(load.average(lambda:load.iterate(filename, oname, hname, filerange, framerange, suffix=suffix), span=avgspan)):
            water = water_module.Molecule()
            # Main part of the program is contained in th Formatter object. (See formats/)
            #logger.debug("Output file format: {0}".format(file_format))
            #formatter = safe_import("format", file_format)
            lattice_info = load.make_lattice_info(oatoms, hatoms, cellmat)
            ice = analice.AnalIce(lattice_info,
                                  comment=comment).analyze_ice(water=water,
                                                               formatter=formatter,
                                                               noise=noise,
                                                               )
            if output is not None:
                sys.stdout = open(output % i, "w")
            sys.stdout.write(ice)
        if stdout is not None:
            # recover stdout
            sys.stdout = stdout


if __name__ == "__main__":
    main()
