desc = {
    "ref": {"0": "Fennell 2005"},
    "usage": "No options available.",
    "brief": 'Hypothetical ice "i".',
}

import genice3.unitcell
import numpy as np
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    i単位胞を定義するクラス。
    """

    def __init__(self):

        waters = np.fromstring(
            """
        0.16666 0.16666 0.0
        0.16666 0.16666 0.47
        0.16666 0.83333 0.0
        0.16666 0.83333 0.47
        0.83333 0.16666 0.0
        0.83333 0.16666 0.47
        0.83333 0.83333 0.0
        0.83333 0.83333 0.47
        0.33333 0.33333 0.235
        0.33333 0.33333 0.705
        0.33333 0.66666 0.235
        0.33333 0.66666 0.705
        0.66666 0.33333 0.235
        0.66666 0.33333 0.705
        0.66666 0.66666 0.235
        0.66666 0.66666 0.705
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "absolute"

        bondlen = 0.4

        density = 0.92

        cell = cellvectors(a=1.0, b=1.0, c=0.94)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            density=density,
            # **kwargs,
        )
