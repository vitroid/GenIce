# Another plan of reactive GenIce3.

import genice2
from genice2.reactive import DependencyCacheMixin, property_depending_on
from genice2.molecules import Molecule
from genice2 import ConfigurationError
from genice2.stage1 import replicate_positions
from genice2.stage2 import grandcell_wrap
from genice2.plugin import safe_import
from genice2.valueparser import put_in_array
from genice2.cell import cellvectors, cellshape
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


@dataclass
class GuestSpec:
    """
    ゲストの情報を表すデータクラス。
    """

    molecule: Molecule
    occupancy: float


def orientations(coord, digraph, cellmat, dopants: dict):
    """
    Does not work when two OHs are colinear
    """

    logger = getLogger()
    # just for a test of pure water
    assert len(coord) == digraph.number_of_nodes(), (
        len(coord),
        digraph.number_of_nodes(),
    )
    if len(dopants):
        logger.info(f"  {dopants} dopants")
    # 通常の氷であればアルゴリズムを高速化できる。

    nnode = len(list(digraph))
    neis = np.zeros([nnode, 2], dtype=int)

    # 仮想ノード用の配列。第0要素は実際には第nnode要素を表す。
    extended_coord = []

    celli = np.linalg.inv(cellmat)
    # v0 = np.zeros([nnode, 3])
    # v1 = np.zeros([nnode, 3])
    for node in digraph:
        if node in dopants:
            h1 = np.array([0.0, 1, 1]) / (2**0.5)
            h2 = np.array([0.0, -1, 1]) / (2**0.5)
            r1 = h1 @ celli
            r2 = h2 @ celli
            # 仮想ノードにさしかえる
            neis[node] = [nnode + len(extended_coord), nnode + len(extended_coord) + 1]
            extended_coord += [coord[node] + r1, coord[node] + r2]
            continue
        succ = list(digraph.successors(node))
        if len(succ) < 2:
            vsucc = (coord[succ] - coord[node]) @ cellmat
            pred = list(digraph.predecessors(node))
            vpred = (coord[pred] - coord[node]) @ cellmat
            vsucc /= np.linalg.norm(vsucc, axis=1)[:, np.newaxis]
            vpred /= np.linalg.norm(vpred, axis=1)[:, np.newaxis]
            if len(vpred) > 2:
                # number of incoming bonds should be <= 2
                vpred = vpred[:2]
            vcomp = assume_tetrahedral_vectors(np.vstack([vpred, vsucc]))
            logger.debug(f"Node {node} vcomp {vcomp} vsucc {vsucc} vpred {vpred}")
            vsucc = np.vstack([vsucc, vcomp])[:2]
            rsucc = vsucc @ celli
            # 仮想ノードにさしかえる
            neis[node] = [nnode + len(extended_coord), nnode + len(extended_coord) + 1]
            extended_coord += [coord[node] + rsucc[0], coord[node] + rsucc[1]]
        else:
            neis[node] = succ

    if len(extended_coord) == 0:
        extended_coord = coord
    else:
        extended_coord = np.vstack([coord, extended_coord])

    # array of donating vectors
    v0 = extended_coord[neis[:, 0]] - coord[:]
    v0 -= np.floor(v0 + 0.5)
    v0 = v0 @ cellmat
    v0 /= np.linalg.norm(v0, axis=1)[:, np.newaxis]
    v1 = extended_coord[neis[:, 1]] - coord[:]
    v1 -= np.floor(v1 + 0.5)
    v1 = v1 @ cellmat
    v1 /= np.linalg.norm(v1, axis=1)[:, np.newaxis]
    # intramolecular axes
    y = v1 - v0
    y /= np.linalg.norm(y, axis=1)[:, np.newaxis]
    z = v0 + v1
    z /= np.linalg.norm(z, axis=1)[:, np.newaxis]
    x = np.cross(y, z, axisa=1, axisb=1)

    rotmatrices = np.zeros([nnode, 3, 3])

    rotmatrices[:, 0, :] = x
    rotmatrices[:, 1, :] = y
    rotmatrices[:, 2, :] = z
    return rotmatrices


def guest_parser(guest_options: list):
    # optionで与えられたguest情報をGuestSpecのリストにする。
    logger = getLogger("guest_parser")
    result = {}
    for guest_option in guest_options:
        cage, guest_specs = guest_option.split("=")
        result[cage] = []
        total_occupancy = 0
        for guest_spec in guest_specs.split("+"):
            if "*" in guest_spec:
                molecule, occupancy = guest_spec.split("*")
                occupancy = float(occupancy)
            else:
                molecule = guest_spec
                occupancy = 1.0
            molecule = safe_import("molecule", molecule).Molecule()
            result[cage].append(GuestSpec(molecule, occupancy))
            total_occupancy += occupancy
        if total_occupancy > 1.0:
            raise ConfigurationError(
                f"Total occupancy of {guest_option} is greater than 1.0"
            )
    logger.debug(f"{result=}")
    return result


def ion_parser(ion_options: list):
    result = {}
    for ion_option in ion_options:
        label, ion_name = ion_option.split("=")
        result[int(label)] = ion_name
    return result


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

            repgraph.add_edge(newi, newj)

    return repgraph


def replicate_fixed_edges(
    repgraph: nx.Graph, fixed: nx.DiGraph, nmol: int
) -> nx.DiGraph:
    # replicate_graphの結果を利用してもっとシンプルに処理したい。
    # repgraph = dg.IceGraph()
    logger = getLogger("replicate_fixed_edges")
    fixedEdges = nx.DiGraph()
    for repi, repj in repgraph.edges():
        i = repi % nmol
        j = repj % nmol
        if fixed.has_edge(i, j):
            fixedEdges.add_edge(repi, repj)
        elif fixed.has_edge(j, i):
            fixedEdges.add_edge(repj, repi)
    for edge in fixedEdges.edges():
        logger.debug(f"* {edge=}")
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
        cell: np.ndarray,
        waters,
        bondlen,
        coord="relative",
        density=None,
        graph=None,
        fixed=None,
        cages=None,
        assess_cages=False,
        shift=(0.0, 0.0, 0.0),
        anions=[],
        cations=[],
    ):

        self.cell = cell

        nmol = len(waters)
        volume = np.linalg.det(self.cell)
        original_density = 18 * nmol / (volume * 1e-21 * 6.022e23)
        self.logger.info(f"{original_density=}")
        if density is not None:
            self.logger.info(f"{density=} specified.")
            scale = (density / original_density) ** (1 / 3)
            self.logger.info(f"{scale=}")
            self.cell /= scale

        self.bondlen = bondlen

        if coord == "absolute":
            self.waters = self.cell.abs2rel(waters)
        else:
            self.waters = waters

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
            self.logger.info("Cages are provided...")
            for pos, label in zip(self.cages[0], self.cages[1]):
                self.logger.info(f"  {label} @ {pos}")
            if assess_cages:
                raise ValueError("Cages cannot be assessed if cages are provided.")
        elif assess_cages:
            self.logger.info("Assessing cages...")
            self.cages = genice2.cage.assess_cages(self.graph, self.waters)
        else:
            self.cages = None
        # cagesは必要ない場合もある。また、格子の情報から推定することもできる。しかし、計算コストは大きいので、明示的に指示して推定すべき。

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
            cages=cages,
            coord=coord,
            bondlen=bondlen,
            # density=density,
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
        "lattice_sites",
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

        a, b, c, A, B, C = cellshape(self.replication_matrix @ self.unitcell.cell)
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

    @property_depending_on("replication_matrix", "unitcell")
    def cell(self):
        # 前後逆かも
        return self.unitcell.cell @ self.replication_matrix

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
    @property_depending_on("graph", "depol_loop", "lattice_sites", "fixedEdges")
    def digraph(self):
        """Makes a directed graph."""
        for edge in self.fixedEdges.edges():
            self.logger.debug(f"+ {edge=}")
        dg = genice_core.ice_graph(
            self.graph,
            vertexPositions=self.lattice_sites,
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
    def lattice_sites(self):
        return replicate_positions(
            self.unitcell.waters, self.replica_vectors, self.replication_matrix
        )

    @property_depending_on("lattice_sites", "unitcell", "replica_vectors")
    def anions(self):
        anion_dict = dict()
        Z = len(self.unitcell.waters)
        for label, ion_name in self.unitcell.anions.items():
            for i in range(len(self.replica_vectors)):
                site = i * Z + label
                anion_dict[site] = ion_name
        return anion_dict

    @property_depending_on("lattice_sites", "unitcell", "replica_vectors")
    def cations(self):
        cation_dict = dict()
        Z = len(self.unitcell.waters)
        for label, ion_name in self.unitcell.cations.items():
            for i in range(len(self.replica_vectors)):
                site = i * Z + label
                cation_dict[site] = ion_name
        return cation_dict

    @property_depending_on("anions", "cations", "lattice_sites")
    def site_occupants(self):
        occupants = ["water"] * len(self.lattice_sites)
        for site, ion_name in self.anions.items():
            occupants[site] = ion_name
        for site, ion_name in self.cations.items():
            occupants[site] = ion_name
        return occupants

    @property_depending_on("graph", "unitcell")
    def fixedEdges(self) -> nx.DiGraph:
        return replicate_fixed_edges(
            self.graph, self.unitcell.fixed, len(self.unitcell.waters)
        )

    @property_depending_on("lattice_sites", "digraph", "unitcell")
    def orientations(self):
        # stages 5
        # the last parameter is self.dopants (if any)
        return orientations(
            self.lattice_sites,
            self.digraph,
            self.unitcell.cell,
            self.anions | self.cations,
        )

    @property_depending_on("unitcell", "replica_vectors", "replication_matrix")
    def cages(self):
        repcagepos = replicate_positions(
            self.unitcell.cages[0], self.replica_vectors, self.replication_matrix
        )
        num_cages_in_unitcell = len(self.unitcell.cages[0])
        repcagetype = [
            self.unitcell.cages[1][i % num_cages_in_unitcell]
            for i in range(len(repcagepos))
        ]
        # unit cellのcagesと同じ構造。
        self.logger.debug(f"{repcagepos=}, {repcagetype=}")
        return repcagepos, repcagetype

    def water_molecules(self, water_model: Molecule):
        mols = {}
        for site in range(len(self.lattice_sites)):
            if "water" == self.site_occupants[site]:
                rel_position = self.lattice_sites[site]
                orientation = self.orientations[site]

                sites = water_model.sites @ orientation + rel_position @ self.cell
                mols[site] = Molecule(
                    name=water_model.name,
                    sites=sites,
                    labels=water_model.labels,
                    is_water=True,
                )
        return mols

    def guest_molecules(self, guests: dict):
        randoms = np.random.random(len(self.cages[0]))

        mols = []
        for pos, label, probability in zip(self.cages[0], self.cages[1], randoms):
            accum = 0.0
            if label in guests:
                for guest_spec in guests[label]:
                    molecule = guest_spec.molecule
                    occupancy = guest_spec.occupancy
                    accum += occupancy
                    if accum > probability:
                        mols.append(
                            Molecule(
                                name=molecule.name,
                                sites=molecule.sites + pos @ self.cell,
                                labels=molecule.labels,
                                is_water=molecule.is_water,
                            )
                        )
                        break
        return mols

    def substitutional_ions(self):
        # 将来は分子イオン(H3O+など)を置く可能性もあることに留意。
        ions = {}
        # ならべかえはここではしない。formatterにまかせる。
        for label, name in self.anions.items():
            ions[label] = Molecule(
                name=name,
                sites=[self.lattice_sites[label] @ self.cell],
                labels=[
                    name,
                ],
                is_water=False,
            )
        for label, name in self.cations.items():
            ions[label] = Molecule(
                name=name,
                sites=[self.lattice_sites[label] @ self.cell],
                labels=[
                    name,
                ],
                is_water=False,
            )
        return ions

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

    def to_gro(self, waters: dict, guests: dict, ions: dict):
        """Output in Gromacs format.

        parametersには、hooksで指定したstageの結果が含まれる。
        """
        logger = getLogger()

        cellmat = self.cell

        if not (cellmat[0, 1] == 0 and cellmat[0, 2] == 0 and cellmat[1, 2] == 0):
            logger.info(
                "  The specified reshaping matrix does not obey the requirements for Gromacs' unit cell convention."
            )
            a = np.linalg.norm(cellmat[0])
            b = np.linalg.norm(cellmat[1])
            c = np.linalg.norm(cellmat[2])
            ea = cellmat[0] / a
            eb = cellmat[1] / b
            ec = cellmat[2] / c
            A = np.degrees(np.arccos(eb @ ec))
            B = np.degrees(np.arccos(ec @ ea))
            C = np.degrees(np.arccos(ea @ eb))
            rotmat = np.linalg.inv(cellmat) @ cellvectors(a, b, c, A, B, C)
            logger.info("  The reshape matrix is reoriented.")
        else:
            rotmat = np.eye(3)

        atoms = []
        residue_count = 1
        for water in waters.values():
            for name, position in zip(water.labels, water.sites):
                atoms.append([residue_count, water.name, name, position])
            residue_count += 1
        for guest in guests:
            for name, position in zip(guest.labels, guest.sites):
                atoms.append([residue_count, guest.name, name, position])
            residue_count += 1
        for ion in ions.values():
            for name, position in zip(ion.labels, ion.sites):
                atoms.append([residue_count, ion.name, name, position])
            residue_count += 1

        logger.info("  Total number of atoms: {0}".format(len(atoms)))
        if len(atoms) > 99999:
            logger.warn(
                "  Fixed-digit format of Gromacs cannot deal with atoms more than 99999. Residue number and atom number are set appropriately."
            )
        s = ""
        s += "Generated by GenIce https://github.com/vitroid/GenIce \n"
        s += "{0}\n".format(len(atoms))
        molorder = 0
        for i, atom in enumerate(atoms):
            resno, resname, atomname, position = atom
            position = position @ rotmat
            if resno == 0:
                molorder += 1
            if len(atoms) > 99999:
                s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                    9999, resname, atomname, 9999, position[0], position[1], position[2]
                )
            else:
                s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                    molorder,
                    resname,
                    atomname,
                    i + 1,
                    position[0],
                    position[1],
                    position[2],
                )
        cellmat = cellmat @ rotmat
        if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
            s += "    {0:.8f} {1:.8f} {2:.8f}\n".format(
                cellmat[0, 0], cellmat[1, 1], cellmat[2, 2]
            )
        else:
            s += "    {0:.8f} {1:.8f} {2:.8f} {3:.8f} {4:.8f} {5:.8f} {6:.8f} {7:.8f} {8:.8f}\n".format(
                cellmat[0, 0],
                cellmat[1, 1],
                cellmat[2, 2],
                cellmat[0, 1],
                cellmat[0, 2],
                cellmat[1, 0],
                cellmat[1, 2],
                cellmat[2, 0],
                cellmat[2, 1],
            )
        return s


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
@click.option(
    "--guest", "-g", multiple=True, help="Guest descriptors, e.g. A12=me, A14=co2, etc."
)
@click.option(
    "--anion",
    "-a",
    multiple=True,
    help="Specify a monatomic anion that replaces a water molecule., e.g. -a 1=Cl, -a 35=Br, etc., where the number is the label of a water molecule in the unit cell.",
)
@click.option(
    "--cation",
    "-c",
    multiple=True,
    help="Specify a monatomic cation that replaces a water molecule., e.g. -c 1=Na, -c 35=K, etc., where the number is the label of a water molecule in the unit cell.",
)
@click.option(
    "--density", "--dens", type=float, default=None, help="Density of the ice"
)
def main(
    debug,
    shift,
    depol_loop,
    replication_matrix,
    replication_factors,
    assess_cages,
    guest,
    anion,
    cation,
    density,
):
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()
    # shift = tuple(shift)
    if replication_matrix is None:
        replication_matrix = np.diag(replication_factors)
    else:
        replication_matrix = np.array(replication_matrix)
    guest_info = guest_parser(guest)
    anion_info = ion_parser(anion)
    cation_info = ion_parser(cation)
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
    genice.unitcell = A15(
        shift=shift,
        assess_cages=assess_cages,
        anions=anion_info,
        cations=cation_info,
        density=density,
    )

    # # stage 2
    # print(genice.graph)
    # # stage 4
    # print(genice.digraph)
    # # genice.unitcell = ice1h(shift=shift)
    # # stage 5
    # # print(genice.orientations)
    # # stage 6 and 7
    # # オプションを指定できるものはmethodとし、指定しないものはpropertyとする。
    # print(genice.molecules(types=[MoleculeType.WATER]))
    waters = genice.water_molecules(water_model=FourSiteWater())
    # あとは
    # guest: ネットワークには影響しないので、moleculesで書き加えればいい。
    guests = genice.guest_molecules(guests=guest_info)
    # spot_guest: その前に、replicate後のcageの場所を何らかの方法でユーザーに知らせないといけない。まああとまわしでいいだろう。
    # anion/cation: fixedEdgesに影響する。
    ions = genice.substitutional_ions()
    # spot_dopants: これも要るだろうね。
    # group

    # 試しにgroファイルにしてみるか。
    with open("genice3.gro", "w") as f:
        f.write(genice.to_gro(waters, guests, ions))


if __name__ == "__main__":
    main()
