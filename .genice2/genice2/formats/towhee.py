# coding: utf-8
"""
towhee file format (.coods?) [EXPERIMENTAL]
"""

import numpy as np
from logging import getLogger
from genice2 import rigid
from genice2.decorators import timeit, banner
import genice2.formats
from genice2.molecules import serialize


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in the Towhee format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {7: self.Hook7}

    @timeit
    @banner
    def Hook7(self, ice):
        "Output in TowHee format."
        logger = getLogger()
        cellmat = ice.repcell.mat

        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)

        s = ""
        for i, atom in enumerate(atoms):
            molorder, resname, atomname, position, order = atom
            s += "{0:20.15f}{1:20.15f}{2:20.15f}   {3:6}{4:3}{5:6}{6:6}\n".format(
                position[0] *
                10,
                position[1] *
                10,
                position[2] *
                10,
                atomname,
                resname,
                order +
                1,
                i +
                1)
        print(s, end="")
        s = "  Cell Dimensions:"
        logger.info(s)
        s = "    a: {0} {1} {2}".format(cellmat[0, 0] * 10,
                                        cellmat[0, 1] * 10,
                                        cellmat[0, 2] * 10)
        logger.info(s)
        s = "    b: {0} {1} {2}".format(cellmat[1, 0] * 10,
                                        cellmat[1, 1] * 10,
                                        cellmat[1, 2] * 10)
        logger.info(s)
        s = "    c: {0} {1} {2}".format(cellmat[2, 0] * 10,
                                        cellmat[2, 1] * 10,
                                        cellmat[2, 2] * 10)
        logger.info(s)
