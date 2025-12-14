#!/usr/bin/python
# coding: utf-8


from genice2.lattices.iceM import Lattice as LatticeM

desc = {
    "ref": {"Md": "Mochizuki 2024"},
    "usage": "No options available.",
    "brief": "A hydrogen-disordered counterpart of ice M.",
    "test": ({"options": "--depol=none"},),
}


class Lattice(LatticeM):
    def __init__(self):
        super().__init__()
        del self.fixed



# ============================================================================
# New genice3.unitcell implementation (TODO: implement manually)
# ============================================================================

desc = {'ref': {'Md': 'Mochizuki 2024'}, 'usage': 'No options available.', 'brief': 'A hydrogen-disordered counterpart of ice M.', 'test': ({'options': '--depol=none'},)}

import genice3.unitcell
import numpy as np
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    iceMd単位胞を定義するクラス。

    NOTE: This unitcell is not yet implemented.
    Please contact the maintainer or implement it manually.
    """

    def __init__(self, **kwargs):
        raise NotImplementedError(
            f"{self.__class__.__name__} is not yet implemented. "
            "This unitcell requires manual implementation. "
            "Please contact the maintainer or implement it manually. "
            f"Reason: watersが定義されていないため, cellが定義されていないため"
        )