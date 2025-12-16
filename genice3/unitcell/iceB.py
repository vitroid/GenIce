#!/usr/bin/python

import numpy as np
from genice3.util import (
    atomdic,
    symmetry_operators,
    waters_and_pairs,
    density_in_g_cm3,
)
import genice3.unitcell
import networkx as nx
from cif2ice import cellvectors

desc = {
    "ref": {"B": "Baez 1998"},
    "usage": "No options available.",
    "brief": "Hypothetical ice B.",
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self, **kwargs):
        atoms = """
        O1 0.2382 0.7981 0.3422
        O2 0.5000 0.5000 0.0446
        H1 -0.2603 0.9742 0.4613
        H2 0.1446 0.8545 0.2013
        H3 0.4108 0.6151 0.1577
        """

        # P2_1 2_1 2
        symops = """
             x            y            z
            1/2+x        1/2-y         -z
            1/2-x        1/2+y         -z
             -x           -y            z
        """

        # in nm
        a, b, c = 0.7176, 0.4428, 0.5040

        cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        # the unit self.cell is too small to handle; multiply (2,2,2)
        waters, fixed = waters_and_pairs(cell, atomd, sops, rep=(2, 2, 2))

        density = density_in_g_cm3(len(waters), cell)

        coord = "relative"

        super().__init__(
            cell=cell,
            waters=waters,
            density=density,
            coord=coord,
            graph=nx.Graph(fixed),
            fixed=nx.DiGraph(fixed),
        )
