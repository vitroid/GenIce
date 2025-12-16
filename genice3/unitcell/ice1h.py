desc = {
    "ref": {},
    "usage": "No options available.",
    "brief": "Most popular Ice I (hexagonal). "
    + "NOTE: Due to a historical reason, the crystal axes of hexagonal ice are exchanged. "
    + "If you want the basal plane to be Z axis, please use 'one[hh]' instead.",
}

import genice3.unitcell
import numpy as np
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice1h単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        1.328 1.802 3.38
        5.267 4.524 1.109
        6.58 5.442 3.365
        5.267 4.542 5.629
        2.623 0.877 5.644
        2.667 5.488 5.625
        5.241 1.756 1.12
        5.241 1.774 5.64
        1.354 4.588 7.888
        1.354 4.57 3.369
        2.667 5.47 1.105
        2.623 0.858 1.124
        6.537 0.831 3.384
        6.537 0.849 7.903
        6.581 5.461 7.884
        1.328 1.82 7.899
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "absolute"

        bondlen = 3.0

        density = 0.92

        cell = cellvectors(a=7.84813412606925, b=7.37735062301457, c=9.06573834219084)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            density=density,
            **kwargs,
        )
