#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice3.util import cellvectors, atomdic, symmetry_operators, waters_and_pairs
import genice3.unitcell
import networkx as nx

desc = {
    "ref": {"VIII": "Kuhs 1998"},
    "usage": "No options available.",
    "brief": "Ice VIII, a hydrogen-ordered counterpart of ice VII.",
    "test": ({"options": "--depol=none"},),
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):
        atoms = """
        O 0   0.25       0.1071(12)
        D 0   0.4157(15) 0.1935(8)
        """

        # space group: I4_1 /a m d No. 141

        symops = """
             x,            y,            z
             -x,          1/2-y,          z
            3/4-y,        1/4+x,        3/4+z
            3/4+y,        3/4-x,        1/4+z
             -x,            y,            z
              x,          1/2-y,          z
            1/4+y,        1/4+x,        3/4+z
            1/4-y,        3/4-x,        1/4+z

             -x,           -y,           -z
              x,          1/2+y,         -z
            3/4+y,        1/4-x,        3/4-z
            3/4-y,        3/4+x,        1/4-z
              x,           -y,           -z
             -x,          1/2+y,         -z
            1/4-y,        1/4-x,        3/4-z
            1/4+y,        3/4+x,        1/4-z

        """.translate(
            {ord(","): ""}
        )

        # add +1/2, +1/2, +1/2
        lines = ""
        for line in symops.split("\n"):
            cols = line.split()
            if len(cols) == 3:
                line = " ".join([x + "+1/2" for x in cols]) + "\n"
            lines += line

        symops += lines

        a = 4.656 / 10.0  # nm
        b = a
        c = 6.775 / 10.0  # nm
        A = 90
        B = 90
        C = 90

        cell = cellvectors(a, b, c, A, B, C)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        waters, fixed = waters_and_pairs(cell, atomd, sops, rep=(2, 2, 2))

        density = 18 * len(waters) / 6.022e23 / (np.linalg.det(cell) * 1e-21)
        coord = "relative"

        super().__init__(
            cell=cell,
            waters=waters,
            density=density,
            coord=coord,
            graph=nx.Graph(fixed),
            fixed=nx.DiGraph(fixed),
        )
