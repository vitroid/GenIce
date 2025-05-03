"""
GenIce class
"""

import itertools as it
from collections import defaultdict
from logging import getLogger
from typing import Type, Union, List, Dict, Set, Optional, Tuple

import numpy as np
import pairlist as pl
import networkx as nx

# from genice2 import cage

# from genice2 import digraph_unused as dg
from genice2.cell import Cell, rel_wrap, cellshape
from genice2.decorators import banner, timeit, debug_args
from genice2.formats import Format
from genice2.lattices import Lattice
import genice_core

# A virtual monatomic molecule
from genice2.molecules import Molecule, arrange, monatom, one
from genice2.plugin import Group, safe_import
from genice2.valueparser import (
    parse_cages,
    parse_pairs,
    # plugin_option_parser,
    put_in_array,
)

from dataclasses import dataclass, field


def assume_tetrahedral_vectors(v):
    """
    Assume missing vectors at a tetrahedral node.

    Given: known vectors.
    Returns: assumed vectors
    """

    assert len(v) > 0

    if len(v) == 3:
        return [-(v[0] + v[1] + v[2])]

    if len(v) == 2:
        y = v[1] - v[0]
        y /= np.linalg.norm(y)
        z = v[1] + v[0]
        z /= np.linalg.norm(z)
        x = np.cross(y, z)
        v2 = (x * 8.0**0.5 - z) / 3.0
        v3 = (-x * 8.0**0.5 - z) / 3.0
        return [v2, v3]

    if len(v) == 1:
        vr = np.random.rand(3)
        vr /= np.linalg.norm(vr)
        z = v[0] / np.linalg.norm(v[0])
        x = np.cross(z, vr)
        y = np.cross(z, x)
        x1 = -x / 2 + 3.0**0.5 * y / 2
        x2 = -x / 2 - 3.0**0.5 * y / 2
        return [x, x1, x2]

    return []


# def ice_rule(dg: nx.DiGraph, strict=False) -> bool:
#     if strict:
#         for node in dg:
#             if dg.in_degree(node) != 2 or dg.out_degree(node) != 2:
#                 return False
#     else:
#         for node in dg:
#             if dg.in_degree(node) > 2 or dg.out_degree(node) > 2:
#                 return False
#     return True


def orientations(coord, digraph, cell, fixed_nodes: set):
    """
    Does not work when two OHs are colinear
    """

    logger = getLogger()
    # just for a test of pure water
    assert len(coord) == digraph.number_of_nodes(), (
        len(coord),
        digraph.number_of_nodes(),
    )
    logger.info(f"{fixed_nodes} fixed nodes")
    # 通常の氷であればアルゴリズムを高速化できる。

    nnode = len(list(digraph))
    neis = np.zeros([nnode, 2], dtype=int)

    # 仮想ノード用の配列。第0要素は実際には第nnode要素を表す。
    extended_coord = []

    # v0 = np.zeros([nnode, 3])
    # v1 = np.zeros([nnode, 3])
    for node in digraph:
        if node in fixed_nodes:
            h1 = np.array([0.0, 1, 1]) / (2**0.5)
            h2 = np.array([0.0, -1, 1]) / (2**0.5)
            r1 = cell.abs2rel(h1)
            r2 = cell.abs2rel(h2)
            # 仮想ノードにさしかえる
            neis[node] = [nnode + len(extended_coord), nnode + len(extended_coord) + 1]
            extended_coord += [coord[node] + r1, coord[node] + r2]
            continue
        succ = list(digraph.successors(node))
        if len(succ) < 2:
            vsucc = cell.rel2abs(rel_wrap(coord[succ] - coord[node]))
            pred = list(digraph.predecessors(node))
            vpred = cell.rel2abs(rel_wrap(coord[pred] - coord[node]))
            vsucc /= np.linalg.norm(vsucc, axis=1)[:, np.newaxis]
            vpred /= np.linalg.norm(vpred, axis=1)[:, np.newaxis]
            if len(vpred) > 2:
                # number of incoming bonds should be <= 2
                vpred = vpred[:2]
            vcomp = assume_tetrahedral_vectors(np.vstack([vpred, vsucc]))
            logger.debug(f"Node {node} vcomp {vcomp} vsucc {vsucc} vpred {vpred}")
            vsucc = np.vstack([vsucc, vcomp])[:2]
            rsucc = cell.abs2rel(vsucc)
            # 仮想ノードにさしかえる
            neis[node] = [nnode + len(extended_coord), nnode + len(extended_coord) + 1]
            extended_coord += [coord[node] + rsucc[0], coord[node] + rsucc[1]]
        else:
            if len(succ) > 2:
                logger.info(succ)
            neis[node] = succ

    if len(extended_coord) == 0:
        extended_coord = coord
    else:
        extended_coord = np.vstack([coord, extended_coord])

    # array of donating vectors
    v0 = extended_coord[neis[:, 0]] - coord[:]
    v0 -= np.floor(v0 + 0.5)
    v0 = v0 @ cell.mat
    v0 /= np.linalg.norm(v0, axis=1)[:, np.newaxis]
    v1 = extended_coord[neis[:, 1]] - coord[:]
    v1 -= np.floor(v1 + 0.5)
    v1 = v1 @ cell.mat
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


def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99

    if pairs is None:
        iter = it.combinations(coord, 2)
    else:
        iter = [(coord[i], coord[j]) for i, j in pairs]

    for c1, c2 in iter:
        r = cell.rel2abs(rel_wrap(c1 - c2))
        rr = r @ r

        if rr < dmin:
            dmin = rr

    return dmin**0.5


def replicate_groups(groups, waters, cagepos, rep):
    """
    This is not that easy.
    """
    logger = getLogger()
    # Storage for replicated groups
    newgroups = defaultdict(dict)

    nmol = len(waters)

    for root, cages in groups.items():
        # Position of root (water) (fractional)
        root_pos = waters[root]

        for cage, group_name in cages.items():
            # Position of the cage (fractional)
            cage_pos = cagepos[cage]
            # Relative position of the cage
            delta = rel_wrap(cage_pos - root_pos)
            # (Image) cell that the cage resides
            gcell = np.floor(root_pos + delta)

            for i, xyz in enumerate(rep):
                newroot = i * nmol + root
                replica_vector_b = xyz + gcell
                fractional_b = replica_vector_b @ np.linalg.inv(grand_cellmat)
                fractional_b -= np.floor(fractional_b)
                replica_vector_b = fractional_b @ grand_cellmat

                # test
                d = replica_vector_b - np.floor(replica_vector_b + 0.5)
                assert d @ d < 1e-20

                replica_vector_b = replica_vector_b.astype(int)

                b = replica_vector_labels[tuple(replica_vector_b)]
                newcage = cage + ncage * b
                newgroups[newroot][newcage] = group_name
                # たぶんこれでいいと思うのだが、この部分は今は使っていないので検証不能。

            # for x in range(rep[0]):
            #     for y in range(rep[1]):
            #         for z in range(rep[2]):
            #             r = np.array((x, y, z))
            #             # label of the root (water) in the replica
            #             newroot = root + len(waters) * (z + rep[2] * (y + rep[1] * x))
            #             # replicated cell in which the cage resides.
            #             # modulo by positive number is always positive.
            #             cr = (r + gcell) % rep
            #             newcage = cage + len(cagepos) * (
            #                 cr[2] + rep[2] * (cr[1] + rep[1] * cr[0])
            #             )
            #             newcage = int(newcage)
            #             newgroups[newroot][newcage] = group_name
            #             # logger.info(("root",newroot,"newcage", newcage))
    return newgroups


def replicate_labeldict(labels, nmol, replica_vectors):
    newlabels = dict()

    for j in labels:
        for i, _ in enumerate(replica_vectors):
            newj = nmol * i + j
            newlabels[newj] = labels[j]

    return newlabels


def replicate_graph(
    graph1,
    cell1frac_coords,
    replica_vectors,
    fixed: list,
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
      fixed (list): グラフ内の固定エッジを表すリストです。
      replica_vector_labels: レプリカベクトル座標のタプルを一意のラベルにマッピングする辞書です。このラベルは、複製されたグラフ内の各レプリカ ベクトルを識別するために使用されます。
      reshape: 単位胞を積みかさねて拡大された結晶構造を作る、積み重ね方を表す行列。

    Returns:
      関数 `replicate_graph` は 2 つの値、`repgraph` と `fixedEdges` を返します。
    """
    # repgraph = dg.IceGraph()
    logger = getLogger()
    repgraph = nx.Graph()
    fixed = nx.DiGraph(fixed)
    fixedEdges = nx.DiGraph()
    nmol = cell1frac_coords.shape[0]

    # 正の行列式の値(倍率)。整数。
    det = np.linalg.det(reshape)
    if det < 0:
        det = -det
    det = np.floor(det + 0.5).astype(int)
    # 逆行列に行列式をかけたもの。整数行列。
    invdet = np.floor(np.linalg.inv(reshape) * det + 0.5).astype(int)

    for i, j in graph1.edges(data=False):
        if fixed.has_edge(i, j):
            edge_fixed = True
        elif fixed.has_edge(j, i):
            edge_fixed = True
            i, j = j, i
        else:
            edge_fixed = False
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

            if edge_fixed:  # fixed pair
                fixedEdges.add_edge(newi, newj)
            repgraph.add_edge(newi, newj)

    return repgraph, fixedEdges


def replicate_positions(positions1, replica_vectors, grand_cellmat):
    # レプリカ単位胞の数だけ、水分子位置を複製する。
    inv = np.linalg.inv(grand_cellmat)
    reppositions = []
    for replica_vector in replica_vectors:
        # レプリカ単位胞内での分子の位置
        replica = positions1 + replica_vector
        # 大セル内の位置に換算する
        replica = replica @ inv
        replica -= np.floor(replica)
        # 束ねて
        reppositions.append(replica)
    # 一つの配列にまとめ
    reppositions = np.vstack(reppositions)
    return reppositions


def neighbor_cages_of_dopants(dopants, waters, cagepos, cell):
    """
    Just shows the environments of the dopants
    """
    # logger = getLogger()
    dnei = defaultdict(set)

    for site, name in dopants.items():
        org = waters[site]

        for i, pos in enumerate(cagepos):
            # Displacement (relative)
            a = cell.rel2abs(rel_wrap(pos - org))
            sqdistance = a @ a

            if sqdistance < 0.57**2:
                dnei[site].add(i)
                # logger.info((i,cagepos[i]))

    return dnei


def grandcell_wrap(
    cell1frac_vec: np.ndarray, reshape: np.ndarray, invdet: np.ndarray, det: int
):
    """
    単位胞の小数座標で表された点(単位胞内とは限らない)を、拡大単位胞内におさめる。
    例えば、拡大単位胞が(-1,4), (6,1)なら、(5.1,5.1)と(6.1,1.1)は(0.1,0.1)になる。
    """
    frac = cell1frac_vec @ invdet
    frac = frac % det
    return frac @ reshape / det


@dataclass
class Stage1Output:
    """Stage1の出力データ"""

    reppositions: np.ndarray  # 複製された分子位置
    repcell: Cell  # 複製されたセル
    repcagetype: Optional[List[str]]  # 複製されたケージタイプ
    repcagepos: Optional[np.ndarray]  # 複製されたケージ位置
    cagetypes: Optional[Dict[str, Set[int]]]  # ケージタイプの集合


@dataclass
class Stage2Output:
    """Stage2の出力データ"""

    graph: nx.Graph  # ネットワークトポロジー
    dopants: Dict[int, str]  # ドーパント情報
    fixed_edges: nx.DiGraph  # 固定エッジ
    fixed_nodes: Set[int]  # ドーパントの原子のインデックスの集合
    groups: Dict[int, Dict[int, str]]  # グループ情報
    filled_cages: Set[int]  # 埋められたケージ


@dataclass
class Stage34EOutput:
    """Stage34Eの出力データ"""

    digraph: nx.DiGraph  # 有向グラフ


@dataclass
class Stage5Output:
    """Stage5の出力データ"""

    rotmatrices: np.ndarray  # 回転行列


@dataclass
class Stage6Output:
    """Stage6の出力データ"""

    universe: List[np.ndarray]  # 原子位置


@dataclass
class Stage7Output:
    """Stage7の出力データ"""

    universe: List[np.ndarray]  # 原子位置（ゲスト分子を含む）


class Stage1:
    """単位胞の複製を行うステージ"""

    def __init__(
        self,
        waters1: np.ndarray,
        replica_vectors: np.ndarray,
        reshape_matrix: np.ndarray,
        cell1: Cell,
        cagepos1: Optional[np.ndarray] = None,
        cagetype1: Optional[List[str]] = None,
    ):
        self.waters1 = waters1
        self.replica_vectors = replica_vectors
        self.reshape_matrix = reshape_matrix
        self.cell1 = cell1
        self.cagepos1 = cagepos1
        self.cagetype1 = cagetype1

    @timeit
    @banner
    def execute(self) -> Stage1Output:
        """Stage 1: Replicates the unit cell."""
        reppositions = replicate_positions(
            self.waters1, self.replica_vectors, self.reshape_matrix
        )
        repcell = Cell(self.reshape_matrix @ self.cell1.mat)

        if self.cagepos1 is not None:
            repcagepos = replicate_positions(
                self.cagepos1, self.replica_vectors, self.reshape_matrix
            )
            repcagetype = [
                self.cagetype1[i % len(self.cagetype1)] for i in range(len(repcagepos))
            ]
            cagetypes = defaultdict(set)
            for i, typ in enumerate(repcagetype):
                cagetypes[typ].add(i)
        else:
            repcagepos = None
            repcagetype = None
            cagetypes = None

        return Stage1Output(
            reppositions=reppositions,
            repcell=repcell,
            repcagetype=repcagetype,
            repcagepos=repcagepos,
            cagetypes=cagetypes,
        )


class Stage2:
    """ランダムグラフの生成と複製を行うステージ"""

    def __init__(
        self,
        graph1: nx.Graph,
        waters1: np.ndarray,
        replica_vectors: np.ndarray,
        fixed1: List[Tuple[int, int]],
        replica_vector_labels: Dict[Tuple[int, int, int], int],
        reshape_matrix: np.ndarray,
        anions: Dict[int, str],
        cations: Dict[int, str],
        groups1: Dict[int, Dict[int, str]],
        cagepos1: np.ndarray,
    ):
        self.graph1 = graph1
        self.waters1 = waters1
        self.replica_vectors = replica_vectors
        self.fixed1 = fixed1
        self.replica_vector_labels = replica_vector_labels
        self.reshape_matrix = reshape_matrix
        self.anions = anions
        self.cations = cations
        self.dopants1 = anions | cations
        self.groups1 = groups1
        self.cagepos1 = cagepos1

    @timeit
    @banner
    def execute(self) -> Stage2Output:
        """Stage 2: Makes a random graph and replicates it."""
        logger = getLogger()

        # # Some edges are directed when ions are doped.
        # if self.doping_hook_function is not None:
        #     self.doping_hook_function(self)  # may be defined in the plugin

        # Replicate the dopants in the unit cell
        dopants = replicate_labeldict(
            self.dopants1, len(self.waters1), self.replica_vectors
        )

        groups = replicate_groups(
            self.groups1, self.waters1, self.cagepos1, self.replica_vectors
        )

        # self.groups_info(self.groups)
        filled_cages = set()
        for _, cages in groups.items():
            filled_cages |= set(cages)

        graph, fixed_edges = replicate_graph(
            self.graph1,
            self.waters1,
            self.replica_vectors,
            self.fixed1,
            self.replica_vector_labels,
            self.reshape_matrix,
        )
        # ドーパントの原子のインデックスの集合
        fixed_nodes = set(dopants)

        for site, _ in self.anions.items():
            for nei in self.graph1[site]:
                fixed_edges.add_edge(nei, site)
        for site, _ in self.cations.items():
            for nei in self.graph1[site]:
                fixed_edges.add_edge(site, nei)

        # dopants = replicate_labeldict(
        #     self.dopants1, len(self.waters1), self.replica_vectors
        # )
        # groups = replicate_groups(
        #     self.groups1, self.waters1, self.cagepos1, self.replica_vectors
        # )
        # filled_cages = set()
        # for _, cages in groups.items():
        #     filled_cages |= set(cages)
        logger.info(f"graph: {self.graph1} {graph}")
        # return Stage2Output(dopants, groups, filled_cages, graph, fixedEdges)
        return Stage2Output(
            graph=graph,
            dopants=dopants,
            fixed_edges=fixed_edges,
            fixed_nodes=fixed_nodes,
            groups=groups,
            filled_cages=filled_cages,
        )


class Stage34E:
    """有向氷グラフの生成を行うステージ"""

    def __init__(
        self, graph: nx.Graph, reppositions: np.ndarray, fixedEdges: nx.DiGraph
    ):
        self.graph = graph
        self.reppositions = reppositions
        self.fixedEdges = fixedEdges

    @timeit
    @banner
    def execute(self, depol: str = "none") -> Stage34EOutput:
        """Stage 3: Makes a directed graph."""
        iter = 0 if depol == "none" else 1000
        digraph = genice_core.ice_graph(
            self.graph.to_undirected(),
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=iter,
            fixedEdges=self.fixedEdges,
        )
        return Stage34EOutput(digraph=digraph)


class Stage5:
    """剛体分子の配向を準備するステージ"""

    def __init__(
        self,
        reppositions: np.ndarray,
        digraph: nx.DiGraph,
        repcell: Cell,
        repimmutables: Set[int],
    ):
        self.reppositions = reppositions
        self.digraph = digraph
        self.repcell = repcell
        self.repimmutables = repimmutables

    @timeit
    @banner
    def execute(self) -> Stage5Output:
        """Stage 5: Prepare orientations for rigid molecules."""
        rotmatrices = orientations(
            self.reppositions, self.digraph, self.repcell, self.repimmutables
        )
        return Stage5Output(rotmatrices=rotmatrices)


class Stage6:
    """水分子と置換分子の原子配置を行うステージ"""

    def __init__(
        self,
        reppositions: np.ndarray,
        repcell: Cell,
        rotmatrices: np.ndarray,
        water: Molecule,
        dopants: Dict[int, str],
    ):
        self.reppositions = reppositions
        self.repcell = repcell
        self.rotmatrices = rotmatrices
        self.water = water
        self.dopants = dopants

    @timeit
    @banner
    def execute(self) -> Stage6Output:
        """Stage 6: Arrange atoms of water and replacements"""
        universe = []
        universe.append(
            arrange(
                self.reppositions,
                self.repcell,
                self.rotmatrices,
                self.water,
                immutables=set(self.dopants),
            )
        )
        return Stage6Output(universe=universe)


class Stage7:
    """ゲスト分子の配置を行うステージ"""

    def __init__(
        self,
        universe: List[np.ndarray],
        repcagepos: np.ndarray,
        repcagetype: List[str],
        cagetypes: Dict[str, Set[int]],
        filled_cages: Set[int],
        guests: Dict[str, Dict[str, float]],
        spot_guests: Dict[int, str],
        spot_groups: Dict[int, str],
        groups: Dict[int, Dict[int, str]],
        dopants: Dict[int, str],
        reppositions: np.ndarray,
        rotmatrices: np.ndarray,
        repcell: Cell,
    ):
        self.universe = universe
        self.repcagepos = repcagepos
        self.repcagetype = repcagetype
        self.cagetypes = cagetypes
        self.filled_cages = filled_cages
        self.guests = guests
        self.spot_guests = spot_guests
        self.spot_groups = spot_groups
        self.groups = groups
        self.dopants = dopants
        self.reppositions = reppositions
        self.rotmatrices = rotmatrices
        self.repcell = repcell

    @timeit
    @banner
    def execute(self) -> Stage7Output:
        """Stage 7: Place guest molecules."""
        # ゲスト分子の配置処理
        # ... (既存のStage7の処理を実装)
        return Stage7Output(self.universe)


@dataclass
class GenIceConfig:
    """GenIceの設定を保持するデータクラス"""

    signature: str = ""
    density: float = 0
    rep: Optional[Tuple[int, int, int]] = None
    reshape: np.ndarray = field(default_factory=lambda: np.eye(3, dtype=int))
    cations: Dict[int, str] = field(default_factory=dict)
    anions: Dict[int, str] = field(default_factory=dict)
    spot_guests: Dict[int, str] = field(default_factory=dict)
    spot_groups: Dict[int, str] = field(default_factory=dict)
    asis: bool = False
    shift: Tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass
class GenIceState:
    """GenIceの状態を保持するデータクラス"""

    waters1: np.ndarray
    cell1: Cell
    pairs1: Optional[List[Tuple[int, int]]] = None
    bondlen: Optional[float] = None
    cagepos1: Optional[np.ndarray] = None
    cagetype1: Optional[List[str]] = None
    fixed1: List[Tuple[int, int]] = field(default_factory=list)
    dopants1: Set[int] = field(default_factory=set)
    # groups1: Dict[int, Dict[int, str]] = field(
    #     default_factory=lambda: defaultdict(dict)
    # )
    filled_cages: Set[int] = field(default_factory=set)
    graph1: nx.Graph = field(default_factory=nx.Graph)


class GenIce:
    """氷の結晶構造を生成するクラス"""

    def __init__(
        self,
        lat: Type[Lattice],
        config: Optional[GenIceConfig] = None,
    ):
        logger = getLogger()
        self.config = config or GenIceConfig()
        self.state = self._initialize_state(lat)
        self._setup_replica_vectors()
        self._setup_documentation(lat)
        self._setup_waters(lat)
        self._setup_bond_length(lat)
        self._setup_pairs(lat)
        self._setup_density(lat)
        self._setup_cages(lat)
        self._setup_fixed_bonds(lat)
        # self._setup_dopants(lat)
        self._setup_graph(lat)
        self._setup_groups(lat)

    def _initialize_state(self, lat: Type[Lattice]) -> GenIceState:
        """初期状態を設定"""
        return GenIceState(
            waters1=np.array([]),
            cell1=Cell(lat.cell),
        )

    def _setup_replica_vectors(self):
        """レプリカベクトルの設定"""
        logger = getLogger()
        if self.config.rep is not None:
            logger.warning("rep for GenIce() is deprecated. Use reshape instead.")
            self.replica_vectors = np.array(
                [
                    (x, y, z)
                    for x in range(self.config.rep[0])
                    for y in range(self.config.rep[1])
                    for z in range(self.config.rep[2])
                ]
            )
            self.reshape_matrix = np.diag(self.config.rep)
        else:
            self._setup_reshape_matrix()

    def _setup_reshape_matrix(self):
        """reshape行列の設定"""
        logger = getLogger()
        logger.info("  Reshaping the unit cell.")
        self.reshape_matrix = self.config.reshape

        i, j, k = np.array(self.reshape_matrix)
        logger.info(f"    i:{i}")
        logger.info(f"    j:{j}")
        logger.info(f"    k:{k}")

        a, b, c, A, B, C = cellshape(self.reshape_matrix @ self.state.cell1.mat)
        logger.info("  Reshaped cell:")
        logger.info(f"    a,b,c = {a}, {b}, {c}")
        logger.info(f"    A,B,C = {A}, {B}, {C}")

        self.replica_vectors = self._calculate_replica_vectors(i, j, k)

        # レプリカ単位胞には0から順に番号ラベルがついている。replica_vector_labelsは単位胞の位置をラベルに変換する
        self.replica_vector_labels = {
            tuple(xyz): i for i, xyz in enumerate(self.replica_vectors)
        }

    def _calculate_replica_vectors(
        self, i: np.ndarray, j: np.ndarray, k: np.ndarray
    ) -> np.ndarray:
        """レプリカベクトルの計算"""
        logger = getLogger()
        corners = np.array(
            [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
        )

        mins = np.min(corners, axis=0)
        maxs = np.max(corners, axis=0)
        logger.info(f"mins: {mins}, maxs: {maxs}")

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

    def _setup_documentation(self, lat: Type[Lattice]):
        """ドキュメントの設定"""
        logger = getLogger()
        try:
            self.doc = lat.__doc__.splitlines()
        except BaseException:
            self.doc = []

        if len(self.config.signature) > 0:
            self.doc.append("")
            self.doc.append(self.config.signature)

        for line in self.doc:
            logger.info("  " + line)

    def _setup_waters(self, lat: Type[Lattice]):
        """水分子の設定"""
        logger = getLogger()
        self.state.waters1 = put_in_array(lat.waters)
        logger.debug(f"Waters: {len(self.state.waters1)}")
        self.state.waters1 = self.state.waters1.reshape((-1, 3))

        if lat.coord == "absolute":
            self.state.waters1 = self.state.cell1.abs2rel(self.state.waters1)

        self.state.waters1 = np.array(self.state.waters1) + np.array(self.config.shift)
        self.state.waters1 -= np.floor(self.state.waters1)

    def _setup_pairs(self, lat: Type[Lattice]):
        """分子間結合の設定"""
        logger = getLogger()

        if "pairs" in lat.__dict__ and lat.pairs is not None:
            self.state.pairs1 = parse_pairs(lat.pairs)
        else:
            logger.info("  Pairs are not given explicitly.")
            logger.info("  Estimating the bonds according to the pair distances.")

            logger.debug(f"Bondlen: {self.state.bondlen}")
            # make bonded pairs according to the pair distance.
            # make before replicating them.
            self.state.pairs1 = [
                (i, j)
                for i, j in pl.pairs_iter(
                    self.state.waters1,
                    self.state.bondlen,
                    self.state.cell1.mat,
                    distance=False,
                )
            ]

    def _setup_bond_length(self, lat: Type[Lattice]):
        """結合長の設定"""
        logger = getLogger()
        nmol = self.state.waters1.shape[0]
        volume = self.state.cell1.volume()

        try:
            self.state.bondlen = lat.bondlen
            logger.info(f"Bond length (specified): {self.state.bondlen}")
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            rc = (volume / nmol) ** (1 / 3) * 1.5
            p = pl.pairs_iter(
                self.state.waters1,
                maxdist=rc,
                cell=self.state.cell1.mat,
                distance=False,
            )
            self.state.bondlen = 1.1 * shortest_distance(
                self.state.waters1, self.state.cell1, pairs=p
            )
            logger.info(f"Bond length (estim.): {self.state.bondlen}")

    def _setup_density(self, lat: Type[Lattice]):
        """密度の設定"""
        logger = getLogger()
        nmol = self.state.waters1.shape[0]
        volume = self.state.cell1.volume()
        mass = 18  # water
        NB = 6.022e23
        density0 = mass * nmol / (NB * volume * 1e-21)

        if self.config.density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                logger.info(
                    "Density is not specified. Assume the density from lattice."
                )
                dmin = shortest_distance(self.state.waters1, self.state.cell1)
                logger.info(
                    f"Closest pair distance: {dmin} (should be around 0.276 nm)"
                )
                self.density = density0 / (0.276 / dmin) ** 3
        else:
            self.density = self.config.density

        logger.info(f"Target Density: {self.density}")
        logger.info(f"Original Density: {density0}")

        ratio = (density0 / self.density) ** (1.0 / 3.0)
        self.state.cell1.scale(ratio)

        if self.state.bondlen is not None:
            self.state.bondlen *= ratio
        logger.info(f"Bond length (scaled, nm): {self.state.bondlen}")

    def _setup_cages(self, lat: Type[Lattice]):
        """ケージの設定"""
        logger = getLogger()
        if "cages" in lat.__dict__:
            self.state.cagepos1, self.state.cagetype1 = parse_cages(lat.cages)
            logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            self.state.cagepos1, self.state.cagetype1 = (
                np.array(lat.cagepos),
                lat.cagetype,
            )

        if self.state.cagepos1 is not None:
            self.state.cagepos1 = np.array(self.state.cagepos1) + np.array(
                self.config.shift
            )
            self.state.cagepos1 -= np.floor(self.state.cagepos1)

    def _setup_fixed_bonds(self, lat: Type[Lattice]):
        """固定結合の設定"""
        logger = getLogger()
        if "fixed" in lat.__dict__ and lat.fixed is not None:
            self.state.fixed1 = parse_pairs(lat.fixed)
            logger.info("Orientations of some edges are fixed.")

        if self.config.asis and len(self.state.fixed1) == 0:
            self.state.fixed1 = self.state.pairs1

    # def __setup_dopants(self, lat: Type[Lattice]):
    #     """Hook関数の設定"""
    #     if "dopeIonsToUnitCell" in lat.__dict__:
    #         self.doping_hook_function = lat.dopeIonsToUnitCell
    #     else:
    #         self.doping_hook_function = None

    def _setup_graph(self, lat: Type[Lattice]):
        """グラフの設定"""
        logger = getLogger()
        logger.info(f"Generating the graph...")

        self.state.graph1 = nx.Graph(self.state.pairs1)
        logger.info(f"graph: {self.state.graph1}")

    def _setup_groups(self, lat: Type[Lattice]):
        """グループの設定"""
        logger = getLogger()
        logger.info("Generating the groups...")
        self.state.groups1 = defaultdict(dict)

    def generate_ice(
        self,
        formatter: Type[Format],
        water: Union[Type[Molecule], None] = None,
        guests={},
        depol="strict",
        noise=0.0,
        assess_cages=False,
    ):
        """氷の結晶構造を生成する"""
        logger = getLogger()

        hooks = formatter.hooks()
        max_stage = max(0, *hooks.keys())
        parameters = dict()
        # Stage1: 水分子の複製
        stage1 = Stage1(
            self.state.waters1,
            self.replica_vectors,
            self.reshape_matrix,
            self.state.cell1,
            self.state.cagepos1,
            self.state.cagetype1,
        )
        stage1_output = stage1.execute()

        if 1 in hooks:
            parameters[1] = stage1_output
            if hooks[1] is not None:
                hooks[1](self.state, parameters)
            if max_stage == 1:
                return formatter.dump()

        # Stage2: ランダムグラフの生成と複製
        stage2 = Stage2(
            self.state.graph1,
            self.state.waters1,
            self.replica_vectors,
            self.state.fixed1,
            self.replica_vector_labels,
            self.reshape_matrix,
            self.config.anions,
            self.config.cations,
            self.state.groups1,
            self.state.cagepos1,
        )
        stage2_output = stage2.execute()

        if 2 in hooks:
            parameters[2] = stage2_output
            if hooks[2] is not None:
                hooks[2](self.state, parameters)
            if max_stage == 2:
                return formatter.dump()

        # Stage34E: 有向氷グラフの生成
        stage34e = Stage34E(
            stage2_output.graph, stage1_output.reppositions, stage2_output.fixed_edges
        )
        stage34e_output = stage34e.execute(depol)

        if 3 in hooks:
            parameters[3] = stage34e_output
            if hooks[3] is not None:
                hooks[3](self.state, parameters)
            if max_stage == 3:
                return formatter.dump()

        if 4 in hooks:
            parameters[4] = stage34e_output
            if hooks[4] is not None:
                hooks[4](self.state, parameters)
            if max_stage == 4:
                return formatter.dump()

        # Stage5: 剛体分子の配向を準備
        stage5 = Stage5(
            stage1_output.reppositions,
            stage34e_output.digraph,
            stage1_output.repcell,
            stage2_output.fixed_nodes,
        )
        stage5_output = stage5.execute()

        if 5 in hooks:
            parameters[5] = stage5_output
            if hooks[5] is not None:
                hooks[5](self.state, parameters)
            if max_stage == 5:
                return formatter.dump()

        # Stage6: 水分子と置換分子の原子配置
        stage6 = Stage6(
            stage1_output.reppositions,
            stage1_output.repcell,
            stage5_output.rotmatrices,
            water,
            stage2_output.dopants,
        )
        stage6_output = stage6.execute()

        if 6 in hooks:
            parameters[6] = stage6_output
            if hooks[6] is not None:
                hooks[6](self.state, parameters)
            if max_stage == 6:
                return formatter.dump()

        # Stage7: ゲスト分子の配置
        stage7 = Stage7(
            stage6_output.universe,
            stage1_output.repcagepos,
            stage1_output.repcagetype,
            stage1_output.cagetypes,
            stage2_output.filled_cages,
            guests,
            self.config.spot_guests,
            self.config.spot_groups,
            stage2_output.groups,
            stage2_output.dopants,
            stage1_output.reppositions,
            stage5_output.rotmatrices,
            stage1_output.repcell,
        )
        stage7_output = stage7.execute()

        if 7 in hooks:
            parameters[7] = stage7_output
            if hooks[7] is not None:
                hooks[7](self.state, parameters)

        return formatter.dump()


def dopants_info(dopants=None, waters=None, cagepos=None, cell=None):
    logger = getLogger()
    if dopants is None:
        dopants = self.dopants

    if waters is None:
        waters = self.waters

    if cagepos is None:
        cagepos = self.cagepos

    if cell is None:
        cell = self.cell

    dopants_neighbors = neighbor_cages_of_dopants(dopants, waters, cagepos, cell)

    for dopant, cages in dopants_neighbors.items():
        logger.info(f"    Cages adjacent to dopant {dopant}: {cages}")

    return dopants_neighbors


def groups_info(groups):
    logger = getLogger()
    for root, cages in groups.items():
        for cage, group in cages.items():
            logger.info(f"    Group {group} of dopant {root} in cage {cage}")


def guests_info(cagetypes, molecules):
    logger = getLogger()
    for cagetype, cageid in cagetypes.items():
        logger.info(f"    Guests in cage type {cagetype}:")

        for molec, cages in molecules.items():
            cages = set(cages)
            cages &= cageid

            if len(cages):
                logger.info(f"      {molec} * {len(cages)} @ {cages}")
