#!/usr/bin/python
# coding: utf-8

import numpy as np
from genice3.util import (
    cellvectors,
    atomdic,
    symmetry_operators,
    waters_and_pairs,
    density_in_g_cm3,
)
import genice3.unitcell
import networkx as nx

desc = {
    "ref": {"A": "Baez 1998"},
    "usage": "No options available.",
    "brief": "Hypothetical ice A.",
}


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self):
        atoms = """
        O1 0.0995 0.1901 0.2835
        O2 0.3431 0.3469 0.7598
        O3 0.4055 0.0558 0.5013
        H1 -0.0149 0.1524 0.3658
        H2 0.0921 0.3365 0.2570
        H3 0.4441 0.3484 0.8630
        H4 0.3713 0.2326 0.6746
        H5 0.3045 0.1078 0.4119
        H6 0.3847 -0.0897 0.5194
        """

        # P4_1, No. 76
        symops = """
              x            y            z
             -x           -y          1/2+z
             -y            x          1/4+z
              y           -x          3/4+z
        """

        # in nm
        a, b, c = 0.6733, 0.6733, 0.7164

        cell = cellvectors(a, b, c)

        # helper routines to make from CIF-like data
        atomd = atomdic(atoms)
        sops = symmetry_operators(symops)
        waters, fixed = waters_and_pairs(cell, atomd, sops)

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
