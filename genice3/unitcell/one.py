# coding: utf-8
"""
Generate a hydrogen-disordered ice I with stacking disorder.

Usage:
  genice2 one[hcchchcc]            Specify layer types with "c" or "h".
  genice2 one[hh]                  Pure hexagonal ice one. (Stacking cycle of 2)
  genice2 one[ccc]                 Pure cubic ice one. (Stacking cycle of 3)
"""


import genice3.unitcell
import numpy as np
from cif2ice import cellvectors
from logging import getLogger


desc = {
    "ref": {},
    "usage": __doc__,
    "brief": "Ice I w/ stacking disorder.",
    "test": (
        {
            "args": "ccchchc",  # argument for the plugin itself
        },
        {
            "args": {"layers": "ccchchc"},  # argument for the plugin
        },
    ),
}


lat = [
    [[0, 0], [2, 0], [1, 3], [3, 3]],
    [[0, 2], [2, 2], [1, 5], [3, 5]],
    [[0, 4], [2, 4], [1, 1], [3, 1]],
]


class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self, **kwargs):
        logger = getLogger()

        layers = ""
        if len(kwargs) == 1:
            if "layers" in kwargs:
                layers = kwargs["layers"]
            # chのみで構成された文字列なら
            else:
                for k, v in kwargs.items():
                    if v is True and isinstance(k, str) and all(c in "CcHh" for c in k):
                        layers = k
                        break
        assert layers != "", "Stacking pattern must be specified."

        layer = 0
        height = 0
        dir = 1
        L = []
        for ch in layers:
            for x, y in lat[layer]:
                L.append([x, y, height])
            layer = (layer + dir + 3) % 3
            height += 1
            for x, y in lat[layer]:
                L.append([x, y, height])
            height += 3
            assert ch in "CcHh"
            if ch in "Hh":
                # hexagonal = alternative
                dir = -dir
                # cubic = progressive
        assert layer == 0 and dir == 1, "Incompatible number of layers."
        assert len(L) > 0, "Stacking pattern must be specified."
        waters = np.array(L) / np.array([4.0, 6.0, height])
        coord = "relative"
        LHB = 0.276
        bondlen = 0.3
        y = LHB * (8**0.5 / 3) * 3
        x = y * 2 / 3**0.5
        z = LHB * height / 3
        cell = cellvectors(x, y, z)
        density = 0.92

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            density=density,
            # **kwargs,
        )
