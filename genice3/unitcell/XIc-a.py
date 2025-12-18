"""
Data sources

"""

desc = {
    "ref": {"XIc": "Geiger 2014"},
    "usage": "No options available.",
    "brief": "A candidate for the proton-ordered counterpart of ice Ic. The structure 'a' in Figure 1.",
}

import genice3.unitcell
import numpy as np
import networkx as nx
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    XIc-a単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

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

        # bondlen = 1.9

        density = 0.92

        cell = cellvectors(a=4.0, b=4.0, c=4.0)

        fixed_pairs = [
            (0, 4),
            (0, 5),
            (1, 4),
            (1, 5),
            (2, 6),
            (2, 7),
            (3, 6),
            (3, 7),
            (4, 2),
            (4, 3),
            (5, 2),
            (5, 3),
            (6, 0),
            (6, 1),
            (7, 0),
            (7, 1),
        ]

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            coord=coord,
            # bondlen=bondlen,
            density=density,
            fixed=nx.DiGraph(fixed_pairs),
            graph=nx.Graph(fixed_pairs),
            **kwargs,
        )
