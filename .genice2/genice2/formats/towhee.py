# coding: utf-8
"""
towhee file format (.coods?) [EXPERIMENTAL]
"""

from logging import getLogger
from io import TextIOWrapper
import sys

from genice2.decorators import timeit, banner
import genice2.formats
from genice2.molecules import serialize
from genice2.genice import GenIce


class Format(genice2.formats.Format):
    """
    The atomic positions of the molecules are output in the Towhee format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output in TowHee format."
        logger = getLogger()
        cellmat = genice.cell_matrix()

        atoms = []
        for pos in genice.full_atomic_positions():
            atoms += serialize(pos)

        s = ""
        for i, atom in enumerate(atoms):
            molorder, resname, atomname, position, order = atom
            s += "{0:20.15f}{1:20.15f}{2:20.15f}   {3:6}{4:3}{5:6}{6:6}\n".format(
                position[0] * 10,
                position[1] * 10,
                position[2] * 10,
                atomname,
                resname,
                order + 1,
                i + 1,
            )
        print(s, end="")
        s = "  Cell Dimensions:"
        logger.info(s)
        s = "    a: {0} {1} {2}".format(
            cellmat[0, 0] * 10, cellmat[0, 1] * 10, cellmat[0, 2] * 10
        )
        logger.info(s)
        s = "    b: {0} {1} {2}".format(
            cellmat[1, 0] * 10, cellmat[1, 1] * 10, cellmat[1, 2] * 10
        )
        logger.info(s)
        s = "    c: {0} {1} {2}".format(
            cellmat[2, 0] * 10, cellmat[2, 1] * 10, cellmat[2, 2] * 10
        )
        logger.info(s)
        file.write(s)
