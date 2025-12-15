desc = {
    "ref": {"VI": "Petrenko 1999"},
    "usage": "No options available.",
    "brief": "Conventional high-pressure ice VI.",
}

import genice3.unitcell
import numpy as np
import networkx as nx
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice6単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        pairs_str = """
        1 0
        3 2
        4 2
        4 3
        5 0
        5 1
        6 0
        6 1
        6 5
        7 2
        7 3
        7 4
        8 0
        8 1
        8 5
        8 6
        9 2
        9 3
        9 4
        9 7
        """.split(
            "\n"
        )
        pairs = []
        for line in pairs_str:
            cols = line.split()
            if len(cols) == 2:
                pairs.append((int(cols[0]), int(cols[1])))

        graph = nx.Graph(pairs)

        waters = np.fromstring(
            """
        0.7504    0.2502    0.7510
        0.4707    0.2502    0.3666
        0.9710    0.7509    0.6348
        0.5298    0.7509    0.6348
        0.2501    0.4710    0.8674
        0.0295    0.2502    0.3666
        0.7504    0.5301    0.1341
        0.2501    0.7509    0.2503
        0.7504    0.9716    0.1341
        0.2501    0.0295    0.8674
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        # bondlen = 3.0

        density = 1.373

        cell = cellvectors(a=6.181, b=6.181, c=5.698)

        super().__init__(
            cell=cell,
            waters=waters,
            graph=graph,
            coord=coord,
            # bondlen=bondlen,
            density=density,
            **kwargs,
        )
