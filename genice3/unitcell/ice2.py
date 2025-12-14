#!/usr/bin/python
"""
Usage: genice3 2
"""

desc = {
    "ref": {"2atom": "Kamb 2003", "2cell": "Kamb 1964", "C1": "Londono 1988"},
    "usage": "\nUsage: genice3 2\n",
    "brief": "Hydrogen-ordered ice II.",
}

import genice3.unitcell
import numpy as np
from logging import getLogger
from genice3.util import cellvectors, atomdic, symmetry_operators, waters_and_pairs
import networkx as nx


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice2単位胞を定義するクラス。

    NOTE: This unitcell is not yet implemented.
    Please contact the maintainer or implement it manually.
    """

    def __init__(self, **kwargs):
        logger = getLogger()

        # Ref. 2atom
        atoms = """
        O1   0.2716    0.0259   -0.1471
        O2   0.4798    0.7571    0.3389
        D1   0.7284    0.4038    0.4034
        D2   0.1491    0.0412   -0.2023
        D3   0.7420    0.1978    0.3708
        D4   0.4232    0.1954   -0.0164
        """

        # Ref. 2cell
        # space group: R-3
        symops = """
          x,            y,            z
          z,            x,            y
          y,            z,            x

         -x,           -y,           -z
         -z,           -x,           -y
         -y,           -z,           -x
        """.replace(
            ",", " "
        )

        # Ref. 2cell
        a = 7.78 / 10.0  # nm
        A = 113.1

        cell = cellvectors(a, a, a, A, A, A)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        waters, fixed = waters_and_pairs(cell, atomd, sops)

        # set pairs in this way for hydrogen-ordered ices.
        fixed = nx.DiGraph(fixed)
        pairs = nx.Graph(fixed)

        density = 18 * len(waters) / 6.022e23 / (np.linalg.det(cell) * 1e-21)

        coord = "relative"
        super().__init__(
            cell=cell,
            waters=waters,
            graph=pairs,
            density=density,
            coord=coord,
            bondlen=0.3,
        )
