desc = {
    "ref": {
        "sVII": "Jeffrey 1984",
        "CS4": "Kosyakov 1999",
        "SOD": "IZA Database",
        "207_1_4435": "Engel 2018",
        "engel01": "Engel 2018",
    },
    "usage": "No options available.",
    "brief": "Hypothetical zeolitic ice"
}

import genice3.unitcell
import numpy as np
from genice2.cell import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    SOD単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.0 5.85484414822 3.90322943215
        3.90322943215 1.95161471607 0.0
        1.95161471607 0.0 3.90322943215
        5.85484414822 3.90322943215 0.0
        0.0 3.90322943215 1.95161471607
        3.90322943215 0.0 5.85484414822
        0.0 3.90322943215 5.85484414822
        3.90322943215 0.0 1.95161471607
        3.90322943215 5.85484414822 0.0
        0.0 1.95161471607 3.90322943215
        5.85484414822 0.0 3.90322943215
        1.95161471607 3.90322943215 0.0
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "absolute"

        bondlen = 3.0

        # density = 0.75396428378

        cell = cellvectors(a=7.8064588643, b=7.8064588643, c=7.8064588643)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )