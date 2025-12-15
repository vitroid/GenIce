desc = {
    "ref": {"0": "Russo 2014"},
    "usage": "No options available.",
    "brief": 'Metastable ice "0".',
}

import genice3.unitcell
import numpy as np
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice0単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        4.89641699678464 3.88902019840311 8.96065434241106
        3.88902019840311 4.89641699678464 6.72310471665148
        0.960541133340528 3.88902019840311 3.73273465605689
        0 0 5.22791968635418
        1.96793793172206 0.960541133340528 6.72310471665148
        2.92847906506258 2.92847906506258 0
        4.89641699678464 1.96793793172206 3.73273465605689
        1.96793793172206 4.89641699678464 1.4951850302973
        0 0 0
        3.88902019840311 0.960541133340528 1.4951850302973
        0.960541133340528 1.96793793172206 8.96065434241106
        2.92847906506258 2.92847906506258 5.22791968635418
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "absolute"

        bondlen = 3.0

        density = 0.92

        cell = cellvectors(a=5.85695813012517, b=5.85695813012517, c=10.4558393727084)

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            density=density,
            **kwargs,
        )
