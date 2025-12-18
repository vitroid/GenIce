"""
Data sources

"""

desc = {
    "ref": {"Ic": "Vos 1993"},
    "usage": "No options available.",
    "brief": "Cubic type of ice I.",
}

import genice3.unitcell
import numpy as np
import networkx as nx
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ice1c単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        pairs_str = """
        0 4
        0 5
        0 6
        0 7
        1 4
        1 5
        1 6
        1 7
        2 4
        2 5
        2 6
        2 7
        3 4
        3 5
        3 6
        3 7
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
        0 0 0
        0.5 0.5 0
        0.5 0 0.5
        0 0.5 0.5
        0.25 0.25 0.25
        0.75 0.75 0.25
        0.75 0.25 0.75
        0.25 0.75 0.75
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        density = 0.92

        cell = cellvectors(a=4.0, b=4.0, c=4.0)

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            graph=graph,
            coord=coord,
            density=density,
            **kwargs,
        )
