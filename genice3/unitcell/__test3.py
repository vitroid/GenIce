# coding: utf-8
"""
Test for generating water positions from cage positions
"""


import genice2.lattices
from genice2 import FrankKasper
from genice2.valueparser import parse_cages
from genice2.cell import cellvectors


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cages = """
        12 0.5 0.5 0.5
        14 0.5 0.0 -0.25
        14 0.0 0.25 0.5
        14 0.75 0.5 0.0
        14 0.5 0.0 0.25
        12 0.0 0.0 0.0
        14 0.25 0.5 0.0
        14 0.0 -0.25 0.5
        """

        self.cell = cellvectors(a=12.747893943706936,
                                b=12.747893943706936,
                                c=12.747893943706936)
        cagepos, cagetype = parse_cages(self.cages)
        self.lattice_sites = [w for w in FrankKasper.toWater(cagepos, self.cell)]
        self.coord = "relative"
        self.density = FrankKasper.estimate_density(
            self.lattice_sites, self.cell, 2.76)
        self.bondlen = 2.76 * 1.2



# ============================================================================
# New genice3.unitcell implementation (TODO: implement manually)
# ============================================================================

"""
Test for generating water positions from cage positions
"""

import genice3.unitcell
import numpy as np
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    __test3単位胞を定義するクラス。

    NOTE: This unitcell is not yet implemented.
    Please contact the maintainer or implement it manually.
    """

    def __init__(self, **kwargs):
        raise NotImplementedError(
            f"{self.__class__.__name__} is not yet implemented. "
            "This unitcell requires manual implementation. "
            "Please contact the maintainer or implement it manually. "
            f"Reason: lattice_sitesが定義されていないため"
        )