from logging import getLogger
import numpy as np
import networkx as nx
import pairlist as pl
from genice3 import ConfigurationError
from genice3.util import assess_cages, shortest_distance, density_in_g_cm3
from typing import Dict, Any
from genice3.molecule.one import Molecule


def ion_processor(arg: dict) -> Dict[int, Molecule]:
    # keyとvalueを変換するのみ
    result: Dict[int, Molecule] = {}
    for label, molecule in arg.items():
        result[int(label)] = molecule
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
        lattice_sites: np.ndarray,
        bondlen: float = None,  # graphが与えられていればbondlenは不要
        coord: str = "relative",
        density: float = None,
        graph: nx.Graph = None,
        fixed: nx.DiGraph = nx.DiGraph(),
        shift: tuple = (0.0, 0.0, 0.0),
        anion: dict = {},
        cation: dict = {},
        assess_cages: bool = False,
        name: str = "",  # dummy
    ):
        anion = ion_processor(anion)
        cation = ion_processor(cation)
        if type(density) == str:
            density = float(density)

        self.cell = cell
        celli = np.linalg.inv(cell)

        # 格子点の位置をfractional coordinateにする。
        if coord == "absolute":
            self.lattice_sites = lattice_sites @ celli
        else:
            self.lattice_sites = lattice_sites

        self.logger.debug(f"  {shift=}")
        self.lattice_sites += np.array(shift)
        self.lattice_sites -= np.floor(self.lattice_sites)

        nmol = len(lattice_sites)
        volume = np.linalg.det(self.cell)
        original_density = density_in_g_cm3(nmol, self.cell)
        self.logger.info(f"{original_density=}")

        # bondlenとgraphを同時に指定した場合はErrorとする。
        if bondlen is not None and graph is not None:
            raise ValueError("bondlen and graph cannot be specified at the same time.")

        if bondlen is None:
            short = shortest_distance(self.lattice_sites, self.cell)
            bondlen = 1.1 * short

            # densityが指定されていない場合は、ここで推定するが、採用はしない。(cellが正しいと信じる)
            if density is None:
                estimated_density = original_density * (short / 0.276) ** 3
                self.logger.info(
                    f"Neither bond length nor density is specified. Estimated density: {estimated_density}"
                )

        if graph is None:
            self.graph = nx.Graph(
                [
                    (i, j)
                    for i, j in pl.pairs_iter(
                        self.lattice_sites, bondlen, self.cell, distance=False
                    )
                ]
            )
            self.logger.info(
                f"The HB graph is generated from the bond length: {bondlen}"
            )
        else:
            self.graph = graph
        self.logger.debug(f"  {self.graph.size()=} {self.graph.number_of_nodes()=}")

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
        if any(label >= len(self.lattice_sites) for label in anion):
            raise ValueError(
                "Anion labels must be less than the number of water molecules."
            )
        if any(label >= len(self.lattice_sites) for label in cation):
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
            # lattice_sitesは既にshift済みなので、ケージ位置も既にshift済みの座標系で計算される
            self.logger.info("Assessing cages...")
            self._cages = assess_cages(self.graph, self.lattice_sites)

            # ケージ情報をログ出力
            if self._cages is not None:
                for i, (pos, label) in enumerate(
                    zip(self._cages.positions, self._cages.specs)
                ):
                    self.logger.info(f"cage {i}: {label} @ {pos}")

        return self._cages
