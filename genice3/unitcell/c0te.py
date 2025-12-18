#!/usr/bin/python
"""
Usage: genice2 c0te
"""

from genice3.util import (
    atomdic,
    symmetry_operators,
    waters_and_pairs,
    density_in_g_cm3,
)
import genice3.unitcell
from logging import getLogger
import numpy as np
import networkx as nx
from cif2ice import cellvectors

desc = {
    "ref": {"C0": "Teeratchanan 2015"},
    "usage": __doc__,
    "brief": "Filled ice C0 by Teeratchanan (Hydrogen-disordered.) (Positions of guests are supplied.)",
    "test": ({"options": "-r 2 2 2"},),
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):
        logger = getLogger()

        # Ref. 2atom
        atoms = """
        O1 0.2342 0.4721 0.8019
        O2 0.7648 0.5306 0.2941
        Ne1 -0.0647 0.7868 0.7669
        """

        # Ref.
        # space group: P3_2
        symops = """
          x,            y,            z
         -y,x-y,z+2/3
         -x+y,-x,z+1/3
        """.replace(
            ",", " "
        )

        # Ref. 2cell
        a = 6.177 / 10.0  # nm
        c = 6.054 / 10.0  # nm
        C = 120.0

        cell = cellvectors(a, a, c, C=C)

        # helper routines to make from CIF-like data

        atomd = atomdic(atoms)
        # atoms = fullatoms(atomd, symmetry_operators(symops))
        sops = symmetry_operators(symops)

        # cagetype = []
        # cagepos = []
        # for atomname, pos in atoms:
        #     if atomname == "Ne1":
        #         self.cagetype.append(atomname)
        #         self.cagepos.append(pos)

        waters, pairs = waters_and_pairs(cell, atomd, sops)

        density = density_in_g_cm3(len(waters), cell)
        coord = "relative"
        super().__init__(
            cell=cell,
            lattice_sites=waters,
            density=density,
            coord=coord,
            bondlen=0.3,
        )
