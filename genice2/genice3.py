# Another plan of reactive Genice3.

from genice2.reactive import DependencyCacheMixin, property_depending_on
from genice2 import ConfigurationError
from genice2.stage1 import replicate_positions
from genice2.stage2 import grandcell_wrap
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


@dataclass
class GenIce3Config:
    shift: tuple[float, float, float] = (0.0, 0.0, 0.0)


class UnitCell:
    """
    単位胞を定義する基底クラス。

    Reactiveではないので、GenIce3に与えたあとで内容を変更しても、GenIce3の依存関係に影響しない。
    もし内容を変更したいなら、新たなUnitCellオブジェクトを作成して、GenIce3に与える。
    """

    # 単位胞タイプごとに必要なパラメータのセットを定義
    # このセットに含まれるパラメータは既知として扱われる
    REQUIRED_CELL_PARAMS: set[str] = set()

    def __init__(
        self,
        cell: Cell,
        waters,
        density,
        bondlen,
        coord="relative",
        graph=None,
        fixed=None,
        **kwargs,
    ):
        self.cell = cell
        self.density = density
        self.bondlen = bondlen

        if coord == "absolute":
            self.waters = self.cell.abs2rel(waters)
        else:
            self.waters = waters
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


class ice1h(UnitCell):
    # UnitCellオブジェクトの具体例。

    logger = getLogger("ice1h")

    def __init__(self):
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
        )

        # self.pairsが


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

    def __init__(
        self,
        unitcell,
        water_model=None,
        depol_loop=1000,
        reshape_matrix=np.eye(3, dtype=int),
        **kwargs,
    ):
        self.unitcell = unitcell
        self.water = safe_import("molecule", water_model)
        self.depol_loop = depol_loop
        self.kwargs = kwargs
        self.config = GenIce3Config(shift=(0.0, 0.0, 0.0))
        # 初期値。setterを使わない。
        self.reshape_matrix = reshape_matrix

        a, b, c, A, B, C = cellshape(self.reshape_matrix @ self.unitcell.cell.mat)
        self.logger.debug("  Reshaped cell:")
        self.logger.debug(f"    {a=:.4f}, {b=:.4f}, {c=:.4f}")
        self.logger.debug(f"    {A=:.3f}, {B=:.3f}, {C=:.3f}")

    @property
    def unitcell(self):
        return self._unitcell

    @unitcell.setter
    def unitcell(self, unitcell):
        self._unitcell = unitcell
        self.logger.debug(f"  {unitcell=}")
        self.logger.debug(f"  {unitcell.waters=}")
        self.logger.debug(f"  {unitcell.graph=}")
        self.logger.debug(f"  {unitcell.fixed=}")

    @property
    def reshape_matrix(self):
        return self._reshape_matrix

    @reshape_matrix.setter
    def reshape_matrix(self, reshape_matrix):
        self._reshape_matrix = reshape_matrix

        i, j, k = np.array(reshape_matrix)
        self.logger.debug(f"    {i=}")
        self.logger.debug(f"    {j=}")
        self.logger.debug(f"    {k=}")

    @property_depending_on("_reshape_matrix")
    def replica_vectors(self) -> np.ndarray:
        """レプリカベクトルの計算"""
        i, j, k = np.array(self.reshape_matrix)
        corners = np.array(
            [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        )

        mins = np.min(corners, axis=0)
        maxs = np.max(corners, axis=0)
        self.logger.debug(f"  {mins=}, {maxs=}")

        det = abs(np.linalg.det(self.reshape_matrix))
        det = np.floor(det + 0.5).astype(int)
        invdet = np.floor(np.linalg.inv(self.reshape_matrix) * det + 0.5).astype(int)

        vecs = set()
        for a in range(mins[0], maxs[0] + 1):
            for b in range(mins[1], maxs[1] + 1):
                for c in range(mins[2], maxs[2] + 1):
                    abc = np.array([a, b, c])
                    rep = grandcell_wrap(abc, self.reshape_matrix, invdet, det).astype(
                        int
                    )
                    if tuple(rep) not in vecs:
                        vecs.add(tuple(rep))

        vecs = np.array(list(vecs))
        vol = abs(np.linalg.det(self.reshape_matrix))
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
        "reshape_matrix",
    )
    def graph(self):
        g = replicate_graph(
            self.unitcell.graph,
            self.unitcell.waters,
            self.replica_vectors,
            self.replica_vector_labels,
            self.reshape_matrix,
        )
        return g

    @property_depending_on("unitcell", "replica_vectors", "reshape_matrix")
    def reppositions(self):
        return replicate_positions(
            self.unitcell.waters, self.replica_vectors, self.reshape_matrix
        )

    @property_depending_on("graph", "fixed")
    def fixedEdges(self):
        return replicate_fixed_edges(self.graph, self.unitcell.fixed)


# clickを用い、-dオプションでデバッグレベルを指定できるようにする。
@click.command()
@click.option("--debug", "-d", is_flag=True, help="Enable debug mode")
def main(debug):
    basicConfig(level=DEBUG if debug else INFO)
    logger = getLogger()
    logger.debug("Debug mode enabled")
    genice = GenIce3(unitcell=ice1h(), water_model="4site")
    print(genice.digraph)


if __name__ == "__main__":
    main()
