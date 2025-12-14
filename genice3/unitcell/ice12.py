desc = {'ref': {'12(1)': 'Lobban 1998', '12(2)': 'Koza 2000'}, 'usage': 'No options available.', 'brief': 'Metastable high-pressure ice XII.'}

import genice3.unitcell
import numpy as np
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice12単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.75 0.13356 0.1875
        0.75 0.13536 0.6875
        0.86464 0.75 0.3125
        0.86464 0.75 0.8125
        0.13536 0.25 0.3125
        0.13536 0.25 0.8125
        0.5 0 0.375
        0.5 0 0.875
        0.25 0.86464 0.1875
        0.25 0.86464 0.6875
        0 0.5 0.125
        0 0.5 0.625
        0.5 0.5 0.25
        0.5 0.5 0.75
        0 0 0
        0 0 0.5
        0.25 0.63536 0.4375
        0.25 0.63536 0.9375
        0.36464 0.25 0.0625
        0.36464 0.25 0.5625
        0.75 0.36564 0.4375
        0.75 0.36564 0.9375
        0.63536 0.75 0.0625
        0.63536 0.75 0.5625
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 3.0

        # density = 1.4397

        cell = cellvectors(a=8.2816, b=8.2816, c=8.0722)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )