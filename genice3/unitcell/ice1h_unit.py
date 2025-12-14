import genice3.unitcell
import numpy as np
from genice2.cell import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice1h_unit単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.125 0.333 0.125
        0.375 0.833 0.125
        0.125 0.667 0.25
        0.375 0.167 0.25
        0.125 0.667 0.625
        0.375 0.167 0.625
        0.375 0.833 0.75
        0.125 0.333 0.75
        0.625 0.333 0.125
        0.875 0.833 0.125
        0.625 0.667 0.25
        0.875 0.167 0.25
        0.625 0.667 0.625
        0.875 0.167 0.625
        0.875 0.833 0.75
        0.625 0.333 0.75
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 3.0

        # density = 0.92

        cell = cellvectors(a=4.5328691711 * 2, b=7.84813412606925, c=7.37735062301457)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )