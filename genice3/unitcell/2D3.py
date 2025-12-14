desc = {'ref': {}, 'usage': 'No options available.', 'brief': 'Trilayer honeycomb ice.'}

import genice3.unitcell
import numpy as np
import networkx as nx
from genice3.util import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    2D3単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        pairs_str = """
        1 0
        1 3
        1 7
        3 5
        3 6
        4 6
        4 5
        5 7
        6 0
        0 2
        7 2
        4 2
        1 9
        5 13
        6 14
        2 10
        9 8
        9 11
        9 15
        11 13
        11 14
        12 14
        12 13
        13 15
        14 8
        8 10
        15 10
        12 10
        8 16
        11 19
        12 20
        15 23
        17 16
        17 19
        17 23
        19 21
        19 22
        20 22
        20 21
        21 23
        22 16
        16 18
        23 18
        20 18
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
        0.000	0.000	0.
        0.750	0.167	0.05
        0.250	0.167	0.05
        0.750	0.500	0.
        0.250	0.500	0.
        0.500	0.667	0.05
        0.000	0.667	0.05
        0.500	0.000	0.
        0.000	0.000	0.25
        0.750	0.167	0.20
        0.250	0.167	0.20
        0.750	0.500	0.25
        0.250	0.500	0.25
        0.500	0.667	0.20
        0.000	0.667	0.20
        0.500	0.000	0.25
        0.000	0.000	0.40
        0.750	0.167	0.45
        0.250	0.167	0.45
        0.750	0.500	0.40
        0.250	0.500	0.40
        0.500	0.667	0.45
        0.000	0.667	0.45
        0.500	0.000	0.40
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 3

        # density = 0.75

        cell = cellvectors(a=6.928203230275509, b=6.0, c=14.14213562373095)

        super().__init__(
            cell=cell,
            waters=waters,
            graph=graph,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )