#!/usr/bin/python
"""
Usage: genice2 c2te
"""
import genice3.unitcell
from logging import getLogger
import numpy as np
from genice3.util import (
    atomdic,
    symmetry_operators,
    waters_and_pairs,
)
from cif2ice import cellvectors


desc = {
    "ref": {"C2": "Teeratchanan 2015"},
    "usage": __doc__,
    "brief": "Filled ice C2 (cubic ice) by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.)",
}


def pick_atoms(atoms, names, repeat=(1, 1, 1)):
    nrep = np.array(repeat)
    for atomname, fracpos in atoms:
        if atomname in names:
            for x in range(repeat[0]):
                for y in range(repeat[1]):
                    for z in range(repeat[2]):
                        yield atomname, (fracpos + np.array([x, y, z])) / nrep


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):
        logger = getLogger()

        # Ref. C2
        atoms = """
    O1 0.0000 0.0000 0.1937
    Ne1 0.0000 0.0000 0.7066
        """

        # Ref. C2
        # space group: I4_1md No. 109
        # http://img.chem.ucl.ac.uk/sgp/large/109az1.htm
        symops = """
          x,            y,            z
         -x,           -y,            z
         -y,          1/2+x,        1/4+z
          y,          1/2-x,        1/4+z
         -x,            y,            z
          x,           -y,            z
          y,          1/2+x,        1/4+z
         -y,          1/2-x,        1/4+z

          x+1/2,            y+1/2,            z+1/2
         -x+1/2,           -y+1/2,            z+1/2
         -y+1/2,          1/2+x+1/2,        1/4+z+1/2
          y+1/2,          1/2-x+1/2,        1/4+z+1/2
         -x+1/2,            y+1/2,            z+1/2
          x+1/2,           -y+1/2,            z+1/2
          y+1/2,          1/2+x+1/2,        1/4+z+1/2
         -y+1/2,          1/2-x+1/2,        1/4+z+1/2
        """.replace(
            ",", " "
        )

        # Ref. C2
        a = 4.409 / 10.0  # nm
        c = 6.251 / 10.0  # nm

        cell = cellvectors(a, a, c)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        waters, _ = waters_and_pairs(cell, atomd, sops)

        coord = "relative"

        # bondlen = shortest_distance(waters, cell)
        # cell *= 0.276 / bondlen
        density = 0.92

        super().__init__(
            cell=cell,
            waters=waters,
            # bondlen=0.276 * 1.3,
            coord=coord,
            density=density,
        )
