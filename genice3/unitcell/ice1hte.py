#!/usr/bin/python
"""
Usage: genice3 ice1hte
"""
from logging import getLogger
import numpy as np
from genice3.util import (
    cellvectors,
    atomdic,
    symmetry_operators,
    waters_and_pairs,
    fullatoms,
    density_in_g_cm3,
)
import genice3.unitcell
import networkx as nx

desc = {
    "ref": {"C0": "Teeratchanan 2015"},
    "usage": "Usage: genice3 ice1hte\n",
    "brief": "Filled ice Ih by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.)",
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
    def __init__(self, **kwargs):
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
        """.replace(
            ",", " "
        )

        # Ref. Ih
        a = 4.568 / 10.0  # nm
        b = 7.980 / 10.0  # nm
        c = 6.894 / 10.0  # nm

        cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data

        atomd = atomdic(atoms)
        atoms = fullatoms(atomd, symmetry_operators(symops))

        cagetype = []
        cagepos = []
        for name, pos in pick_atoms(atoms, ("Ne1",), repeat=(2, 1, 1)):
            cagetype.append(name)
            cagepos.append(pos)

        waters, pairs = waters_and_pairs(
            cell, atomd, symmetry_operators(symops), rep=(2, 1, 1)
        )

        density = density_in_g_cm3(len(waters), cell)
        coord = "relative"
        super().__init__(
            cell=cell,
            waters=waters,
            density=density,
            coord=coord,
            bondlen=0.3,
        )
