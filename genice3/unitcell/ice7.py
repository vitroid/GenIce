desc = {
    "ref": {},
    "usage": "No options available.",
    "brief": "Conventional high-pressure ice VII.",
}

import genice3.unitcell
import numpy as np
import networkx as nx
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice7単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        pairs_str = """
        0 4
        0 6
        1 4
        1 5
        1 6
        2 6
        3 4
        3 5
        3 6
        3 7
        4 2
        5 0
        5 2
        7 0
        7 1
        7 2
        8 12
        8 14
        9 12
        9 13
        9 14
        10 14
        11 12
        11 13
        11 14
        11 15
        12 10
        13 8
        13 10
        15 8
        15 9
        15 10
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
        0.00 0.00 0.00
        0.50 0.50 0.00
        0.50 0.00 0.50
        0.00 0.50 0.50
        0.25 0.25 0.25
        0.75 0.75 0.25
        0.75 0.25 0.75
        0.25 0.75 0.75
        0.50 0.00 0.00
        0.00 0.50 0.00
        0.00 0.00 0.50
        0.50 0.50 0.50
        0.75 0.25 0.25
        0.25 0.75 0.25
        0.25 0.25 0.75
        0.75 0.75 0.75
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        # bondlen = 3

        density = 1.6

        cell = cellvectors(a=4.0, b=4.0, c=4.0)

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            graph=graph,
            coord=coord,
            # bondlen=bondlen,
            density=density,
            **kwargs,
        )
