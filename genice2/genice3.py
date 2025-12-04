# Another plan of reactive GenIce3.

import genice2
from genice2.reactive import DependencyCacheMixin, property_depending_on
from genice2.molecules import Molecule
from genice2 import ConfigurationError
from genice2.stage1 import replicate_positions
from genice2.stage2 import grandcell_wrap
from genice2.stage5 import orientations
from genice2.plugin import safe_import
from genice2.valueparser import put_in_array
from genice2.cell import Cell, cellvectors, cellshape
import genice_core
import pairlist as pl
import networkx as nx
import numpy as np
from dataclasses import dataclass
from logging import getLogger, DEBUG, INFO, WARNING, ERROR, CRITICAL, basicConfig
from typing import Optional
import click
from enum import Enum
from math import cos, radians, sin
import inspect


# enumのようなものはどう定義する?
class MoleculeType(Enum):
    WATER = "water"
    GUEST = "guest"
    DOPANT = "dopant"
    GROUP = "group"


def replicate_graph(
    graph1,
    cell1frac_coords,
    replica_vectors,
    replica_vector_labels,
    reshape,
):
    """
    関数 `replicate_graph` は、グラフに関連するさまざまな入力を受け取り、指定されたレプリカ ベクトルと形状に基づいてそれを複製し、複製されたグラフと固定エッジを返します。

    2つの座標系がいりみだれているので注意。
    cell1frac: 複製前の単位胞における小数座標
    grandfrac: 複製後の大きな単位胞における小数座標

    Args:
      graph1: 元のグラフを表すグラフ オブジェクトです。
      cell1frac_coords: 元のグラフ内の原子の小数点座標を含む numpy 配列です。
      replica_vectors: レプリカ(拡大単位胞を構成する、もとの単位胞のグリッド)の方向を定義する整数ベクトルのリスト。
      replica_vector_labels: レプリカベクトル座標のタプルを一意のラベルにマッピングする辞書です。このラベルは、複製されたグラフ内の各レプリカ ベクトルを識別するために使用されます。
      reshape: 単位胞を積みかさねて拡大された結晶構造を作る、積み重ね方を表す行列。

    Returns:
      関数 `replicate_graph` は `repgraph` を返します。

    Stage2のもとの関数と違い、fixedの複製は行いません。
    """
    # repgraph = dg.IceGraph()
    repgraph = nx.Graph()
    nmol = cell1frac_coords.shape[0]

    # 正の行列式の値(倍率)。整数。
    det = np.linalg.det(reshape)
    if det < 0:
        det = -det
    det = np.floor(det + 0.5).astype(int)
    # 逆行列に行列式をかけたもの。整数行列。
    invdet = np.floor(np.linalg.inv(reshape) * det + 0.5).astype(int)

    for i, j in graph1.edges(data=False):
        # positions in the original small cell
        cell1_delta = cell1frac_coords[j] - cell1frac_coords[i]
        cell1_delta = np.floor(cell1_delta + 0.5).astype(int)

        for a, cell1frac_a in enumerate(replica_vectors):
            cell1frac_b = cell1frac_a + cell1_delta
            cell1frac_b = np.floor(
                grandcell_wrap(cell1frac_b, reshape, invdet, det)
            ).astype(int)
            b = replica_vector_labels[tuple(cell1frac_b)]
            newi = nmol * b + i
            newj = nmol * a + j

            repgraph.add_edge(newi, newj, original=(i, j))

    return repgraph


def replicate_fixed_edges(
    repgraph,
    fixed: list,
):
    # replicate_graphの結果を利用してもっとシンプルに処理したい。
    # repgraph = dg.IceGraph()
    fixedEdges = nx.DiGraph()
    for repi, repj, (i, j) in repgraph.edges(data="original"):
        if fixed.has_edge(i, j):
            fixedEdges.add_edge(repi, repj)
        elif fixed.has_edge(j, i):
            fixedEdges.add_edge(repj, repi)
    return fixedEdges


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
        cell: Cell,
        waters,
        density,
        bondlen,
        coord="relative",
        graph=None,
        fixed=None,
        cages=None,
        assess_cages=False,
        **kwargs,
    ):

        self.cell = cell
        self.density = density
        self.bondlen = bondlen

        if coord == "absolute":
            self.waters = self.cell.abs2rel(waters)
        else:
            self.waters = waters

        # shiftはkwargに入っているので、ここで抽出する。
        shift = kwargs.pop("shift", (0.0, 0.0, 0.0))
        self.logger.debug(f"  {shift=}")
        self.waters += np.array(shift)
        self.waters -= np.floor(self.waters)
        # もしkwargsに"pairs"がなければ、watersなどから再計算する。
        if graph is None:
            self.graph = nx.Graph(
                [
                    (i, j)
                    for i, j in pl.pairs_iter(
                        self.waters, self.bondlen, self.cell.mat, distance=False
                    )
                ]
            )
        else:
            self.graph = graph

        self.fixed = nx.DiGraph(fixed)
        if cages:
            self.cages = cages
            if assess_cages:
                raise ValueError("Cages cannot be assessed if cages are provided.")
        elif assess_cages:
            self.cages = genice2.cage.assess_cages(self.graph, self.waters)
            self.logger.info(f"  {self.cages=}")
        else:
            self.cages = None
        # cagesは必要ない場合もある。また、格子の情報から推定することもできる。しかし、計算コストは大きいので、明示的に指示して推定すべき。


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
            cell=Cell(
                desc=cellvectors(
                    a=7.84813412606925, b=7.37735062301457, c=9.06573834219084
                )
            ),
            density=0.92,
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
                cagepos.append(np.array(cols[1:]))
                cagelabel.append(cols[0])
        cages = (cagepos, cagelabel)

        bondlen = 3

        density = 0.6637037332735554

        cell = cellvectors(
            a=12.747893943706936, b=12.747893943706936, c=12.747893943706936
        )
        super().__init__(
            cell=Cell(desc=cell),
            waters=waters,
            graph=nx.Graph(pairs),
            # cages=cages,
            coord=coord,
            bondlen=bondlen,
            density=density,
            **kwargs,
        )


class FourSiteWater(Molecule):
    """
    4サイトモデルの水分子を定義するクラス。
    """

    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.15 / 10
        theta = radians(104.52)

        hy = L1 * sin(theta / 2)
        hz = L1 * cos(theta / 2)
        mz = L2
        sites = np.array(
            [[0.0, 0.0, 0.0], [0.0, hy, hz], [0.0, -hy, hz], [0.0, 0.0, mz]]
        )
        sites -= (sites[1] + sites[2] + sites[3] * 0) / 18
        labels = ["OW", "HW1", "HW2", "MW"]
        name = "ICE"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)


class GenIce3(DependencyCacheMixin):
    # __init__で与える基本的な情報以外を可能な限りreactiveにする。
    # たとえば、digraphが必要になったら、それを生成するために必要な変数をどんどんさかのぼって計算していく。
    # 依存関係を指示し、一旦計算したpropertyは変更がなければcacheして利用する。
    # たとえばこんな感じかな。
    # genice = GenIce(lattice="Ih", water="TIP4P")
    # digraph = genice.digraph
    # geniceの引数も省略可能とするが、不可欠であればエラーとする。

    # Class名でlog表示したい。
    logger = getLogger("GenIce3")

    # ユーザー向けAPIとして公開するpropertyのリスト
    # このリストに含まれるpropertyのみAPIドキュメントを作成する
    PUBLIC_API_PROPERTIES = [
        "digraph",
        "graph",
        "reppositions",
        "orientations",
        "molecules",
        "unitcell",
        "replication_matrix",
        "depol_loop",
    ]

    def __init__(
        self,
        depol_loop=1000,
        replication_matrix=np.eye(3, dtype=int),
        **kwargs,
    ):
        # Default値が必要なもの
        self.depol_loop = depol_loop
        self.replication_matrix = replication_matrix

        # Default値が不要なもの
        for key in self.list_settable_reactive_properties():
            if key in kwargs:
                setattr(self, key, kwargs.pop(key))
        if kwargs:
            raise ConfigurationError(f"Invalid keyword arguments: {kwargs}.")

    @property
    def depol_loop(self):
        return self._depol_loop

    @depol_loop.setter
    def depol_loop(self, depol_loop):
        self._depol_loop = depol_loop
        self.logger.debug(f"  {depol_loop=}")

    @property
    def unitcell(self):
        if not hasattr(self, "_unitcell") or self._unitcell is None:
            raise ConfigurationError("Unitcell is not set.")
        return self._unitcell

    @unitcell.setter
    def unitcell(self, unitcell):
        self._unitcell = unitcell
        self.logger.debug(f"  {unitcell=}")
        self.logger.debug(f"  {unitcell.waters=}")
        self.logger.debug(f"  {unitcell.graph=}")
        self.logger.debug(f"  {unitcell.fixed=}")

        a, b, c, A, B, C = cellshape(self.replication_matrix @ self.unitcell.cell.mat)
        self.logger.debug("  Reshaped cell:")
        self.logger.debug(f"    {a=:.4f}, {b=:.4f}, {c=:.4f}")
        self.logger.debug(f"    {A=:.3f}, {B=:.3f}, {C=:.3f}")

    @property
    def replication_matrix(self):
        return self._replication_matrix

    @replication_matrix.setter
    def replication_matrix(self, replication_matrix):
        self._replication_matrix = replication_matrix

        i, j, k = np.array(replication_matrix)
        self.logger.debug(f"    {i=}")
        self.logger.debug(f"    {j=}")
        self.logger.debug(f"    {k=}")

    @property_depending_on("replication_matrix")
    def replica_vectors(self) -> np.ndarray:
        """レプリカベクトルの計算"""
        i, j, k = np.array(self.replication_matrix)
        corners = np.array(
            [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        )

        mins = np.min(corners, axis=0)
        maxs = np.max(corners, axis=0)
        self.logger.debug(f"  {mins=}, {maxs=}")

        det = abs(np.linalg.det(self.replication_matrix))
        det = np.floor(det + 0.5).astype(int)
        invdet = np.floor(np.linalg.inv(self.replication_matrix) * det + 0.5).astype(
            int
        )

        vecs = set()
        for a in range(mins[0], maxs[0] + 1):
            for b in range(mins[1], maxs[1] + 1):
                for c in range(mins[2], maxs[2] + 1):
                    abc = np.array([a, b, c])
                    rep = grandcell_wrap(
                        abc, self.replication_matrix, invdet, det
                    ).astype(int)
                    if tuple(rep) not in vecs:
                        vecs.add(tuple(rep))

        vecs = np.array(list(vecs))
        vol = abs(np.linalg.det(self.replication_matrix))
        assert np.allclose(vol, len(vecs)), (vol, vecs)

        return vecs

    @property_depending_on("replica_vectors")
    def replica_vector_labels(self):
        return {tuple(xyz): i for i, xyz in enumerate(self.replica_vectors)}

    # 例えば、digraphをreactiveにするにはどう書けばいい?
    @property_depending_on("graph", "depol_loop", "reppositions", "fixedEdges")
    def digraph(self):
        """Makes a directed graph."""
        dg = genice_core.ice_graph(
            self.graph,
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=self.depol_loop,
            fixedEdges=self.fixedEdges,
        )
        if not dg:
            raise ConfigurationError("Failed to generate a directed graph.")
        return dg

    @property_depending_on(
        "unitcell",
        "replica_vectors",
        "replica_vector_labels",
        "replication_matrix",
    )
    def graph(self):
        g = replicate_graph(
            self.unitcell.graph,
            self.unitcell.waters,
            self.replica_vectors,
            self.replica_vector_labels,
            self.replication_matrix,
        )
        return g

    @property_depending_on("unitcell", "replica_vectors", "replication_matrix")
    def reppositions(self):
        return replicate_positions(
            self.unitcell.waters, self.replica_vectors, self.replication_matrix
        )

    @property_depending_on("graph", "unitcell")
    def fixedEdges(self):
        return replicate_fixed_edges(self.graph, self.unitcell.fixed)

    @property_depending_on("reppositions", "digraph", "unitcell")
    def orientations(self):
        # stages 5
        # the last parameter is self.dopants (if any)
        return orientations(self.reppositions, self.digraph, self.unitcell.cell, set())

    def molecules(self, water_model=None, types=None):
        mols = []
        if MoleculeType.WATER in types:
            if water_model is None:
                raise ConfigurationError("Water model is not set.")
            for rel_position, orientation in zip(self.reppositions, self.orientations):
                sites = (
                    water_model.sites @ orientation
                    + rel_position @ self.unitcell.cell.mat
                )
                mols.append(
                    Molecule(
                        name=water_model.name,
                        sites=sites,
                        labels=water_model.labels,
                        is_water=True,
                    )
                )
        return mols

    @classmethod
    def get_public_api_properties(cls):
        """
        ユーザー向けAPIとして公開されているpropertyのリストを返す。

        Returns:
            list: 公開API property名のリスト
        """
        return cls.PUBLIC_API_PROPERTIES.copy()

    @classmethod
    def list_all_reactive_properties(cls):
        """
        すべてのreactive property（@property_depending_onでデコレートされたもの）を列挙する。

        Returns:
            dict: property名をキー、propertyオブジェクトを値とする辞書
        """
        return {
            name: prop
            for name, prop in inspect.getmembers(
                cls, lambda x: isinstance(x, property_depending_on)
            )
        }

    @classmethod
    def list_public_reactive_properties(cls):
        """
        ユーザー向けAPIとして公開されているreactive propertyのみを列挙する。

        Returns:
            dict: property名をキー、propertyオブジェクトを値とする辞書
        """
        all_reactive = cls.list_all_reactive_properties()
        public_names = set(cls.PUBLIC_API_PROPERTIES)
        return {
            name: prop for name, prop in all_reactive.items() if name in public_names
        }

    # setterを持つreactiveな変数を列挙。
    @classmethod
    def list_settable_reactive_properties(cls):
        return {
            name: prop
            for name, prop in inspect.getmembers(
                cls, lambda x: isinstance(x, property) and x.fset is not None
            )
        }

    @classmethod
    def list_public_settable_reactive_properties(cls):
        return {
            name: prop
            for name, prop in cls.list_settable_reactive_properties().items()
            if name in cls.PUBLIC_API_PROPERTIES
        }


# clickを用い、-dオプションでデバッグレベルを指定できるようにする。
@click.command()
@click.help_option()
@click.option("--debug", "-D", is_flag=True, help="Enable debug mode")
# tupleをどうやってオプションで指定するの? clickのドキュメントを見る。
# click.Tuple(int)を使う。
@click.option(
    "--shift",
    "-S",
    type=click.Tuple([float, float, float]),
    default=(0.0, 0.0, 0.0),
    help="Shift the unit cell",
)
@click.option("--depol_loop", type=int, default=1000, help="Depolarization loop")
# 3x3行列の対角項3つか、9要素全部を指定する。
@click.option(
    "--replication_matrix",
    type=click.Tuple([int, int, int, int, int, int, int, int, int]),
    default=None,
    help="Replication matrix (9 elements)",
)
# 単位胞をx,y,z軸方向に複製する場合の、それぞれの軸方向の個数を表すのに、どんな名称のオプションが良いか。
@click.option(
    "--replication_factors",
    "--rep",
    type=click.Tuple([int, int, int]),
    default=(1, 1, 1),
    help="Replication factors (3 elements)",
)
@click.option("--assess_cages", "-A", is_flag=True, help="Assess cages")
def main(
    debug, shift, depol_loop, replication_matrix, replication_factors, assess_cages
):
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()
    # shift = tuple(shift)
    if replication_matrix is None:
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix)
    logger.debug("Debug mode enabled")
    genice = GenIce3(
        depol_loop=depol_loop,
        replication_matrix=replication_matrix,
    )
    logger.info("Reactive properties:")
    logger.info(f"     All: {genice.list_all_reactive_properties().keys()}")
    logger.info(f"  Public: {genice.list_public_reactive_properties().keys()}")
    logger.info("Settabe reactive properties:")
    logger.info(f"     All: {genice.list_settable_reactive_properties().keys()}")
    logger.info(f"  Public: {genice.list_public_settable_reactive_properties().keys()}")
    # genice.unitcell = Ice1h(shift=shift, assess_cages=assess_cages)
    genice.unitcell = A15(shift=shift, assess_cages=assess_cages)
    # unitcell内のcageを指定したが、まだreplicateしていない。

    # # stage 2
    print(genice.graph)
    # # stage 4
    print(genice.digraph)
    # # genice.unitcell = ice1h(shift=shift)
    # # stage 5
    # # print(genice.orientations)
    # # stage 6 and 7
    # # オプションを指定できるものはmethodとし、指定しないものはpropertyとする。
    # # print(genice.molecules(types=[MoleculeType.WATER]))
    # print(genice.molecules(water_model=FourSiteWater(), types=[MoleculeType.WATER]))
    # あとは
    # guest
    # anion/cation
    # group


if __name__ == "__main__":
    main()
