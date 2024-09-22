#!/usr/bin/python
"""
Usage: genice2 ice1hte
"""
import genice2.lattices
from logging import getLogger
import numpy as np


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {"ref": {"C0": "Teeratchanan 2015"},
        "usage": usage(),
        "brief": "Filled ice Ih by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.)"
        }


def pick_atoms(atoms, names, repeat=(1, 1, 1)):
    nrep = np.array(repeat)
    for atomname, fracpos in atoms:
        if atomname in names:
            for x in range(repeat[0]):
                for y in range(repeat[1]):
                    for z in range(repeat[2]):
                        yield atomname, (fracpos + np.array([x, y, z])) / nrep


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        logger = getLogger()

        # Ref. Ih
        atoms = """
        O1 0.0000 0.6699 0.0488
        O2 0.0000 0.3377 -0.0557
        Ne1 0.0000 0.0013 0.7539
        """

        # Ref. Ih
        # space group: Cmc2_1
        symops = """
          x,            y,            z
         -x,            y,            z
          x,           -y,          1/2+z
         -x,           -y,          1/2+z
          x+1/2,            y+1/2,            z
         -x+1/2,            y+1/2,            z
          x+1/2,           -y+1/2,          1/2+z
         -x+1/2,           -y+1/2,          1/2+z
        """.replace(',', ' ')

        # Ref. Ih
        a = 4.568 / 10.0  # nm
        b = 7.980 / 10.0  # nm
        c = 6.894 / 10.0  # nm

        from genice2.cell import cellvectors
        self.cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        from genice2 import CIF
        atomd = CIF.atomdic(atoms)
        atoms = CIF.fullatoms(atomd, CIF.symmetry_operators(symops))

        self.cagetype = []
        self.cagepos = []
        for name, pos in pick_atoms(atoms, ("Ne1",), repeat=(2, 1, 1)):
            self.cagetype.append(name)
            self.cagepos.append(pos)

        self.waters, self.pairs = CIF.waters_and_pairs(
            self.cell, atomd, CIF.symmetry_operators(symops), rep=(2, 1, 1))

        self.density = 18 * len(self.waters) / 6.022e23 / \
            (np.linalg.det(self.cell) * 1e-21)
        self.coord = "relative"
