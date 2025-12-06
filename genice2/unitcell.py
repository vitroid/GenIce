"""
UnitCellクラスとそのサブクラスを定義するモジュール。
"""

import genice2
from genice2.cell import cellvectors
from genice2 import ConfigurationError
import networkx as nx
import numpy as np
import pairlist as pl
from logging import getLogger


class UnitCell:
    """
    単位胞を定義する基底クラス。

    Reactiveではないので、GenIce3に与えたあとで内容を変更しても、GenIce3の依存関係に影響しない。
    もし内容を変更したいなら、新たなUnitCellオブジェクトを作成して、GenIce3に与える。
    """

    # 単位胞タイプごとに必要なパラメータのセットを定義
    # このセットに含まれるパラメータは既知として扱われる
    REQUIRED_CELL_PARAMS: set[str] = set()

    logger = getLogger()

    def __init__(
        self,
        cell: np.ndarray,
        waters: np.ndarray,
        bondlen: float,
        coord: str = "relative",
        density: float = None,
        graph: nx.Graph = None,
        fixed: nx.DiGraph = nx.DiGraph(),
        cages: tuple = None,
        assess_cages: bool = False,
        shift: tuple = (0.0, 0.0, 0.0),
        anions: dict = {},
        cations: dict = {},
    ):

        self.cell = cell

        # 格子点の位置をfractional coordinateにする。
        if coord == "absolute":
            self.waters = self.cell.abs2rel(waters)
        else:
            self.waters = waters

        self.logger.debug(f"  {shift=}")
        self.waters += np.array(shift)
        self.waters -= np.floor(self.waters)

        # もしkwargsに"pairs"がなければ、watersなどから再計算する。
        # この時はまだ与えられたセルサイズのままで作業する。
        if graph is None:
            self.graph = nx.Graph(
                [
                    (i, j)
                    for i, j in pl.pairs_iter(
                        self.waters, bondlen, self.cell, distance=False
                    )
                ]
            )
        else:
            self.graph = graph

        nmol = len(waters)
        volume = np.linalg.det(self.cell)
        original_density = 18 * nmol / (volume * 1e-21 * 6.022e23)
        self.logger.info(f"{original_density=}")

        # 密度を指定した場合は、セルサイズを調整する。
        if density is not None:
            self.logger.info(f"{density=} specified.")
            scale = (density / original_density) ** (1 / 3)
            self.logger.info(f"{scale=}")
            self.cell /= scale

        self.fixed = fixed

        if cages is not None and assess_cages:
            raise ValueError("Cages cannot be assessed if cages are provided.")

        if cages is not None:
            self.cages = cages
            self.logger.info("Cages are provided...")
        elif assess_cages:
            self.logger.info("Assessing cages...")
            self.cages = genice2.cage.assess_cages(self.graph, self.waters)
        else:
            self.cages = None

        # ケージ位置をshiftして、fractional coordinateにする。
        if self.cages is not None:
            _pos = self.cages[0]
            _pos += np.array(shift)
            _pos -= np.floor(_pos)
            self.cages = (_pos, self.cages[1])
            for pos, label in zip(self.cages[0], self.cages[1]):
                self.logger.info(f"  {label} @ {pos}")

        # anion, cationは単位胞内でのイオンの位置を示すので、番号が単位胞の水分子数未満でなければならない。
        if any(label >= len(self.waters) for label in anions):
            raise ValueError(
                "Anion labels must be less than the number of water molecules."
            )
        if any(label >= len(self.waters) for label in cations):
            raise ValueError(
                "Cation labels must be less than the number of water molecules."
            )
        self.anions = anions
        self.cations = cations
        # ionは水素結合の向きを固定する。
        for label in anions:
            for nei in self.graph[label]:
                if self.fixed.has_edge(label, nei):
                    raise ConfigurationError(f"Impossible to dope an anion at {label}.")
                else:
                    self.fixed.add_edge(nei, label)
        for label in cations:
            for nei in self.graph[label]:
                if self.fixed.has_edge(nei, label):
                    raise ConfigurationError(f"Impossible to dope a cation at {label}.")
                else:
                    self.fixed.add_edge(label, nei)
        for edge in self.fixed.edges():
            self.logger.debug(f"  {edge=}")


class Ice1h(UnitCell):
    # UnitCellオブジェクトの具体例。

    logger = getLogger("ice1h")

    def __init__(self, **kwargs):
        #
        waters = np.fromstring(
            """
        1.328 1.802 3.38
        5.267 4.524 1.109
        6.58 5.442 3.365
        5.267 4.542 5.629
        2.623 0.877 5.644
        2.667 5.488 5.625
        5.241 1.756 1.12
        5.241 1.774 5.64
        1.354 4.588 7.888
        1.354 4.57 3.369
        2.667 5.47 1.105
        2.623 0.858 1.124
        6.537 0.831 3.384
        6.537 0.849 7.903
        6.581 5.461 7.884
        1.328 1.82 7.899
        """,
            sep=" ",
        ).reshape(-1, 3)

        super().__init__(
            cell=cellvectors(
                a=7.84813412606925, b=7.37735062301457, c=9.06573834219084
            ),
            # density=0.92,
            bondlen=3,
            waters=waters,
            coord="absolute",
            **kwargs,
        )


class A15(UnitCell):
    """
    A15単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        pairs_str = """
        17 26
        34 39
        6 10
        33 39
        7 12
        28 35
        32 40
        34 41
        42 30
        43 31
        19 15
        21 44
        3 27
        7 36
        29 4
        23 45
        14 18
        9 18
        25 33
        37 29
        5 19
        8 0
        20 5
        11 1
        27 38
        37 45
        16 3
        25 43
        28 43
        19 27
        21 29
        7 38
        36 16
        32 30
        28 27
        17 38
        30 2
        40 23
        33 20
        26 12
        25 9
        26 11
        22 34
        38 4
        14 28
        6 24
        4 41
        0 21
        7 14
        36 23
        5 17
        25 2
        13 39
        22 35
        14 40
        21 11
        5 2
        32 26
        20 10
        6 4
        6 3
        36 29
        24 39
        16 15
        20 41
        18 44
        23 42
        44 31
        8 1
        15 30
        9 1
        0 10
        13 3
        24 37
        34 45
        0 16
        33 8
        9 32
        31 10
        35 2
        45 1
        18 37
        13 42
        8 42
        11 41
        12 44
        12 15
        13 35
        22 40
        17 22
        24 43
        19 31
        """.split(
            "\n"
        )
        pairs = []
        for line in pairs_str:
            cols = line.split()
            if len(cols) == 2:
                pairs.append((int(cols[0]), int(cols[1])))

        waters = np.fromstring(
            """
        0.625 0.8125 0.5
        0.3125 0.6875 0.3125
        0.3125 0.3125 0.6875
        0.875 0.0 0.6875
        0.6875 0.875 0.0
        0.5 0.375 0.8125
        0.8125 0.8125 0.8125
        0.8125 0.1875 0.1875
        0.375 0.8125 0.5
        0.1875 0.5 0.375
        0.6875 0.6875 0.6875
        0.5 0.625 0.1875
        0.6875 0.3125 0.3125
        0.125 0.0 0.6875
        0.0 0.3125 0.125
        0.625 0.1875 0.5
        0.75 0.0 0.5
        0.5 0.25 0.0
        0.0 0.5 0.25
        0.6875 0.3125 0.6875
        0.5 0.625 0.8125
        0.6875 0.6875 0.3125
        0.3125 0.125 0.0
        0.125 0.0 0.3125
        0.0 0.6875 0.875
        0.1875 0.5 0.625
        0.5 0.375 0.1875
        0.8125 0.1875 0.8125
        0.0 0.3125 0.875
        0.8125 0.8125 0.1875
        0.375 0.1875 0.5
        0.8125 0.5 0.625
        0.3125 0.3125 0.3125
        0.3125 0.6875 0.6875
        0.3125 0.875 0.0
        0.1875 0.1875 0.8125
        0.875 0.0 0.3125
        0.0 0.6875 0.125
        0.6875 0.125 0.0
        0.1875 0.8125 0.8125
        0.1875 0.1875 0.1875
        0.5 0.75 0.0
        0.25 0.0 0.5
        0.0 0.5 0.75
        0.8125 0.5 0.375
        0.1875 0.8125 0.1875
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        cages_str = """
        12 0.5 0.5 0.5
        14 0.5 0.0 -0.25
        14 0.0 0.25 0.5
        14 0.75 0.5 0.0
        14 0.5 0.0 0.25
        12 0.0 0.0 0.0
        14 0.25 0.5 0.0
        14 0.0 -0.25 0.5
        """

        cagepos = []
        cagelabel = []
        for line in cages_str.split("\n"):
            cols = line.split()
            if len(cols) > 0:
                cagepos.append([float(x) for x in cols[1:]])
                cagelabel.append(cols[0])
        cagepos = np.array(cagepos)
        cagepos -= np.floor(cagepos)
        cages = (cagepos, cagelabel)

        bondlen = 3

        # density = 0.6637037332735554

        cell = cellvectors(
            a=1.2747893943706936, b=1.2747893943706936, c=1.2747893943706936
        )
        super().__init__(
            cell=cell,
            waters=waters,
            graph=nx.Graph(pairs),
            # cages=cages,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )
