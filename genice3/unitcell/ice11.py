desc = {
    "ref": {"11#19": "Fan 2010", "11": "Jackson 1997"},
    "usage": "No options available.",
    "brief": "A candidate for an antiferroelectric Ice XI #19.",
}

import genice3.unitcell
import numpy as np
import networkx as nx
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice11単位胞を定義するクラス。
    """

    def __init__(self):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.25 0     0.125
        0.25 0.333 0
        0.25 0     0.5
        0.25 0.333 0.625
        0.75 0     0.125
        0.75 0.333 0
        0.75 0     0.5
        0.75 0.333 0.625
        0    0.5   0.125
        0    0.833 0
        0    0.5   0.5
        0    0.833 0.625
        0.5  0.5   0.125
        0.5  0.833 0
        0.5  0.5   0.5
        0.5  0.833 0.625
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        a = 4.4923 / 10 * 2
        b = 7.7808 / 10
        c = 7.3358 / 10
        cell = cellvectors(a, b, c)

        fixed_pairs = [
            (0, 1),
            (0, 13),
            (1, 3),
            (1, 12),
            (2, 11),
            (2, 0),
            (3, 2),
            (3, 10),
            (4, 5),
            (4, 9),
            (5, 8),
            (5, 7),
            (6, 15),
            (6, 4),
            (7, 6),
            (7, 14),
            (8, 1),
            (8, 10),
            (9, 8),
            (9, 0),
            (10, 11),
            (10, 7),
            (11, 6),
            (11, 9),
            (12, 5),
            (12, 14),
            (13, 12),
            (13, 4),
            (14, 15),
            (14, 3),
            (15, 2),
            (15, 13),
        ]

        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            fixed=nx.DiGraph(fixed_pairs),
            graph=nx.Graph(fixed_pairs),
        )
