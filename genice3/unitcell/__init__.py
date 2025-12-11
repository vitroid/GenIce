from logging import getLogger
import numpy as np
import networkx as nx
import pairlist as pl
from genice3 import ConfigurationError
from genice3.util import assess_cages
from typing import Dict, Any
from genice3.molecule.one import Molecule


def ion_processor(arg: dict) -> Dict[int, Molecule]:
    # keyとvalueを変換するのみ
    result: Dict[int, Molecule] = {}
    for label, molecule in arg.items():
        result[int(label)] = Molecule(name=molecule, label=molecule)
    return result


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
        shift: tuple = (0.0, 0.0, 0.0),
        anion: dict = {},
        cation: dict = {},
        **kwargs,
    ):
        anion = ion_processor(anion)
        cation = ion_processor(cation)
        if type(density) == str:
            density = float(density)

        self.cell = cell
        celli = np.linalg.inv(cell)

        # 格子点の位置をfractional coordinateにする。
        if coord == "absolute":
            self.waters = waters @ celli
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

        # ケージの調査は遅延評価（reactive）にする
        self._cages = None  # 遅延評価用

        # anion, cationは単位胞内でのイオンの位置を示すので、番号が単位胞の水分子数未満でなければならない。
        if any(label >= len(self.waters) for label in anion):
            raise ValueError(
                "Anion labels must be less than the number of water molecules."
            )
        if any(label >= len(self.waters) for label in cation):
            raise ValueError(
                "Cation labels must be less than the number of water molecules."
            )
        self.anions = anion
        self.cations = cation
        # ionは水素結合の向きを固定する。
        for label in anion:
            for nei in self.graph[label]:
                if self.fixed.has_edge(label, nei):
                    raise ConfigurationError(f"Impossible to dope an anion at {label}.")
                else:
                    self.fixed.add_edge(nei, label)
        for label in cation:
            for nei in self.graph[label]:
                if self.fixed.has_edge(nei, label):
                    raise ConfigurationError(f"Impossible to dope a cation at {label}.")
                else:
                    self.fixed.add_edge(label, nei)
        for edge in self.fixed.edges():
            self.logger.debug(f"  {edge=}")

    @property
    def cages(self):
        """
        ケージ位置とタイプを取得する（遅延評価）。

        初回アクセス時にケージの調査が行われる。
        """
        if self._cages is None:
            # ケージの調査が必要になったときに実行
            # watersは既にshift済みなので、ケージ位置も既にshift済みの座標系で計算される
            self.logger.info("Assessing cages...")
            self._cages = assess_cages(self.graph, self.waters)

            # ケージ情報をログ出力
            if self._cages is not None:
                for i, (pos, label) in enumerate(
                    zip(self._cages.positions, self._cages.specs)
                ):
                    self.logger.info(f"cage {i}: {label} @ {pos}")

        return self._cages
