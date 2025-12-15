#!/usr/bin/python

import numpy as np
from genice3.util import cellvectors, atomdic, symmetry_operators, waters_and_pairs
import genice3.unitcell
import networkx as nx

desc = {
    "ref": {"engel17": "Engel 2018", "DDR": "IZA Database"},
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice",
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):
        atoms = """
    O1       0.7281    0.0539    0.0698
    O2       0.1309    0.2618    0.1075
    O3       0.1999    0.3999    0.1718
    O4       0.1187    0.2374    0.2309
    O5       0.2256    0.0000    0.0000
    O6       0.0000    0.0000    0.1956
    O7       0.0000    0.0000    0.1166
        """

        symops = """
+x,+y,+z
2/3+x,1/3+y,1/3+z
1/3+x,2/3+y,2/3+z
-y,+x-y,+z
2/3-y,1/3+x-y,1/3+z
1/3-y,2/3+x-y,2/3+z
-x+y,-x,+z
2/3-x+y,1/3-x,1/3+z
1/3-x+y,2/3-x,2/3+z
-y,-x,+z
2/3-y,1/3-x,1/3+z
1/3-y,2/3-x,2/3+z
-x+y,+y,+z
2/3-x+y,1/3+y,1/3+z
1/3-x+y,2/3+y,2/3+z
+x,+x-y,+z
2/3+x,1/3+x-y,1/3+z
1/3+x,2/3+x-y,2/3+z
-x,-y,-z
2/3-x,1/3-y,1/3-z
1/3-x,2/3-y,2/3-z
+y,-x+y,-z
2/3+y,1/3-x+y,1/3-z
1/3+y,2/3-x+y,2/3-z
+x-y,+x,-z
2/3+x-y,1/3+x,1/3-z
1/3+x-y,2/3+x,2/3-z
+y,+x,-z
2/3+y,1/3+x,1/3-z
1/3+y,2/3+x,2/3-z
+x-y,-y,-z
2/3+x-y,1/3-y,1/3-z
1/3+x-y,2/3-y,2/3-z
-x,-x+y,-z
2/3-x,1/3-x+y,1/3-z
1/3-x,2/3-x+y,2/3-z
        """.replace(
            ",", " "
        )

        a = 13.7950 / 10.0  # nm
        b = 13.7950 / 10.0  # nm
        c = 40.7500 / 10.0  # nm
        A = 90
        B = 90
        C = 120

        cell = cellvectors(a, b, c, A, B, C)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        waters, pairs = waters_and_pairs(cell, atomd, sops)

        density = 18 * len(waters) / 6.022e23 / (np.linalg.det(cell) * 1e-21)

        coord = "relative"

        super().__init__(
            cell=cell,
            waters=waters,
            # graph=nx.Graph(pairs),
            # fixed=nx.DiGraph(fixed),
            bondlen=0.32,
            density=density,
            coord=coord,
        )
