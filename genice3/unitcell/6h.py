desc = {
    "ref": {},
    "usage": "No options available.",
    "brief": "Half lattice of ice VI.",
    "test": ({"options": "-r 2 2 2"},),
}

import genice3.unitcell
import numpy as np
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    6h単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.2200    0.5000    0.3800
        0.7800    0.5000    0.3800
        0.5000    0.2200    0.6200
        0.5000    0.5000    0.0000
        0.5000    0.7800    0.6200
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 2.3681227356441177

        density = 0.6865

        cell = cellvectors(a=4.87672629, b=4.87385128, c=4.49131038)

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            coord=coord,
            bondlen=bondlen,
            density=density,
            **kwargs,
        )
