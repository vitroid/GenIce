# Another plan of reactive GenIce3 (リアクティブプロパティによる実装).

from genice3 import ConfigurationError
from genice3.molecule import Molecule
from genice3.util import (
    replicate_positions,
    grandcell_wrap,
    assume_tetrahedral_vectors,
    CageSpecs,
)
from cif2ice import cellshape, cellvectors
import genice_core
import networkx as nx
import numpy as np
from dataclasses import dataclass
from logging import getLogger, DEBUG, INFO, WARNING, ERROR, CRITICAL, basicConfig
from typing import Optional, Dict, Tuple, List, Any
from enum import Enum
import inspect
import json

from genice3.dependencyengine import DependencyEngine
from genice3.unitcell import UnitCell


class ShowUsageError(Exception):
    """Usage表示を要求する例外

    Args:
        flag_name: フラグ名（例: "?", "help?", "cage?"）
    """

    def __init__(self, flag_name: str, message: str = ""):
        self.flag_name = flag_name
        super().__init__(message or f"Show usage for flag: {flag_name}")


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

    def __repr__(self) -> str:
        return (
            f"GuestSpec(molecule={self.molecule.name!r}, "
            f"occupancy={self.occupancy:.3f})"
        )


@dataclass
class AtomicStructure:
    """
    原子構造データを統合的に保持するデータクラス。
    exporterプラグインが使用するための統一インターフェース。
    """

    waters: Dict[int, Molecule]
    guests: List[Molecule]
    ions: Dict[int, Molecule]
    cell: np.ndarray

    def __repr__(self) -> str:
        return (
            f"AtomicStructure(n_waters={len(self.waters)}, "
            f"n_guests={len(self.guests)}, "
            f"n_ions={len(self.ions)}, "
            f"cell_shape={self.cell.shape})"
        )


def _assume_water_orientations(
    coord: np.ndarray, digraph: nx.DiGraph, cellmat: np.ndarray, dopants: Dict[int, str]
) -> np.ndarray:
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


def _replicate_graph(
    graph1: nx.Graph,
    cell1frac_coords: np.ndarray,
    replica_vectors: np.ndarray,
    replica_vector_labels: Dict[Tuple[int, ...], int],
    reshape: np.ndarray,
) -> nx.Graph:
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


def _replicate_fixed_edges(
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


# ============================================================================
# DependencyEngineタスク関数定義: 依存関係は引数名から自動推論される
# ============================================================================

_genice3_logger = getLogger("GenIce3")


def cell(unitcell: UnitCell, replication_matrix: np.ndarray) -> np.ndarray:
    """拡大単位胞のセル行列"""
    return unitcell.cell @ replication_matrix


def replica_vectors(replication_matrix: np.ndarray) -> np.ndarray:
    """レプリカベクトルを計算する。

    拡大単位胞を構成するために必要な、元の単位胞のグリッド位置を表す
    整数ベクトルのリストを生成します。各ベクトルは、拡大単位胞内での
    単位胞の相対位置を表します。

    Args:
        replication_matrix: 単位胞を複製するための3x3整数行列

    Returns:
        np.ndarray: レプリカベクトルの配列（Nx3配列、Nは拡大単位胞内の単位胞数）
    """
    i, j, k = np.array(replication_matrix)
    corners = np.array(
        [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    )

    mins = np.min(corners, axis=0)
    maxs = np.max(corners, axis=0)
    _genice3_logger.debug(f"  {mins=}, {maxs=}")

    det = abs(np.linalg.det(replication_matrix))
    det = np.floor(det + 0.5).astype(int)
    invdet = np.floor(np.linalg.inv(replication_matrix) * det + 0.5).astype(int)

    vecs = set()
    for a in range(mins[0], maxs[0] + 1):
        for b in range(mins[1], maxs[1] + 1):
            for c in range(mins[2], maxs[2] + 1):
                abc = np.array([a, b, c])
                rep = grandcell_wrap(abc, replication_matrix, invdet, det).astype(int)
                if tuple(rep) not in vecs:
                    vecs.add(tuple(rep))

    vecs = np.array(list(vecs))
    vol = abs(np.linalg.det(replication_matrix))
    assert np.allclose(vol, len(vecs)), (vol, vecs)
    return vecs


def replica_vector_labels(replica_vectors: np.ndarray) -> Dict[Tuple[int, ...], int]:
    """レプリカベクトルラベルの辞書を生成する。

    各レプリカベクトル（整数タプル）を一意のインデックスにマッピングします。
    このマッピングは、拡大単位胞内での単位胞の位置を効率的に管理するために使用されます。

    Args:
        replica_vectors: レプリカベクトルの配列

    Returns:
        Dict[Tuple[int, ...], int]: レプリカベクトル座標タプルからインデックスへのマッピング
    """
    return {tuple(xyz): i for i, xyz in enumerate(replica_vectors)}


def graph(
    unitcell: UnitCell,
    replica_vectors: np.ndarray,
    replica_vector_labels: Dict[Tuple[int, ...], int],
    replication_matrix: np.ndarray,
) -> nx.Graph:
    """拡大単位胞に対応するグラフを生成する。

    基本単位胞のグラフ（水分子間の水素結合ネットワーク）を、
    拡大単位胞全体に複製して統合したグラフを生成します。
    このグラフは、拡大単位胞内のすべての水分子間の接続関係を表します。

    Args:
        unitcell: 基本単位胞オブジェクト
        replica_vectors: レプリカベクトルの配列
        replica_vector_labels: レプリカベクトルラベルの辞書
        replication_matrix: 単位胞を複製するための3x3整数行列

    Returns:
        nx.Graph: 拡大単位胞全体の水素結合ネットワークを表す無向グラフ
    """
    g = _replicate_graph(
        unitcell.graph,
        unitcell.waters,
        replica_vectors,
        replica_vector_labels,
        replication_matrix,
    )
    return g


def lattice_sites(
    unitcell: UnitCell,
    replica_vectors: np.ndarray,
    replication_matrix: np.ndarray,
) -> np.ndarray:
    """格子サイト位置を拡大単位胞全体に複製する。

    基本単位胞内の水分子の座標を、拡大単位胞全体に複製します。
    各水分子の位置は、単位胞の周期境界条件に従って配置されます。

    Args:
        unitcell: 基本単位胞オブジェクト
        replica_vectors: レプリカベクトルの配列
        replication_matrix: 単位胞を複製するための3x3整数行列

    Returns:
        np.ndarray: 拡大単位胞内のすべての水分子の座標（Nx3配列、Nは水分子数）
    """
    return replicate_positions(unitcell.waters, replica_vectors, replication_matrix)


def anions(
    unitcell: UnitCell, replica_vectors: np.ndarray, spot_anions: Dict[int, str]
) -> Dict[int, str]:
    """アニオンイオンの配置情報を拡大単位胞全体に複製する。

    基本単位胞内で定義されたアニオンイオンと、spot_anionsで指定された
    特定位置のアニオンを統合し、拡大単位胞全体でのアニオン配置を返します。

    Args:
        unitcell: 基本単位胞オブジェクト
        replica_vectors: レプリカベクトルの配列
        spot_anions: 特定の格子サイト位置に配置するアニオンの辞書（サイトインデックス -> イオン名）

    Returns:
        Dict[int, str]: 拡大単位胞全体でのアニオン配置（サイトインデックス -> イオン名）
    """
    anion_dict: Dict[int, str] = {}
    Z = len(unitcell.waters)
    for label, ion_name in unitcell.anions.items():
        for i in range(len(replica_vectors)):
            site = i * Z + label
            anion_dict[site] = ion_name
    for label, ion_name in spot_anions.items():
        anion_dict[label] = ion_name
    return anion_dict


def cations(
    unitcell: UnitCell, replica_vectors: np.ndarray, spot_cations: Dict[int, str]
) -> Dict[int, str]:
    """カチオンイオンの配置情報を拡大単位胞全体に複製する。

    基本単位胞内で定義されたカチオンイオンと、spot_cationsで指定された
    特定位置のカチオンを統合し、拡大単位胞全体でのカチオン配置を返します。

    Args:
        unitcell: 基本単位胞オブジェクト
        replica_vectors: レプリカベクトルの配列
        spot_cations: 特定の格子サイト位置に配置するカチオンの辞書（サイトインデックス -> イオン名）

    Returns:
        Dict[int, str]: 拡大単位胞全体でのカチオン配置（サイトインデックス -> イオン名）
    """
    cation_dict: Dict[int, str] = {}
    Z = len(unitcell.waters)
    for label, ion_name in unitcell.cations.items():
        for i in range(len(replica_vectors)):
            site = i * Z + label
            cation_dict[site] = ion_name
    for label, ion_name in spot_cations.items():
        cation_dict[label] = ion_name
    return cation_dict


def site_occupants(
    anions: Dict[int, str], cations: Dict[int, str], lattice_sites: np.ndarray
) -> List[str]:
    """各格子サイトの占有種（水分子またはイオン）のリストを生成する。

    各格子サイトが水分子、アニオン、カチオンのいずれで占有されているかを
    判定し、サイトインデックス順にリストとして返します。

    Args:
        anions: アニオン配置の辞書
        cations: カチオン配置の辞書
        lattice_sites: 格子サイト位置の配列

    Returns:
        List[str]: 各サイトの占有種のリスト（"water"またはイオン名）
    """
    occupants = ["water"] * len(lattice_sites)
    for site, ion_name in anions.items():
        occupants[site] = ion_name
    for site, ion_name in cations.items():
        occupants[site] = ion_name
    return occupants


def fixedEdges(
    graph: nx.Graph,
    unitcell: UnitCell,
    spot_anions: Dict[int, str],
    spot_cations: Dict[int, str],
) -> nx.DiGraph:
    """固定エッジ（水素結合の方向が固定されたエッジ）を拡大単位胞全体に複製する。

    基本単位胞で定義された固定エッジと、spot_anions/spot_cationsで指定された
    イオン位置に基づく固定エッジを統合し、拡大単位胞全体での固定エッジを返します。
    固定エッジは、水素結合の方向が既に決定されている（プロトンが固定されている）
    エッジを表します。

    Args:
        graph: 拡大単位胞全体のグラフ
        unitcell: 基本単位胞オブジェクト
        spot_anions: 特定位置のアニオン配置
        spot_cations: 特定位置のカチオン配置

    Returns:
        nx.DiGraph: 拡大単位胞全体での固定エッジを表す有向グラフ
    """
    dg = _replicate_fixed_edges(graph, unitcell.fixed, len(unitcell.waters))
    for label in spot_anions:
        for nei in graph.neighbors(label):
            if dg.has_edge(label, nei):
                raise ConfigurationError(
                    f"矛盾する辺の固定 at {label}; すでに({label} --> {nei})が固定されています。"
                )
            else:
                dg.add_edge(nei, label)
    for label in spot_cations:
        for nei in graph.neighbors(label):
            if dg.has_edge(nei, label):
                raise ConfigurationError(
                    f"矛盾する辺の固定 at {label}; すでに({nei} --> {label})が固定されています。"
                )
            else:
                dg.add_edge(label, nei)
    return dg


def digraph(
    graph: nx.Graph,
    depol_loop: int,
    lattice_sites: np.ndarray,
    fixedEdges: nx.DiGraph,
) -> nx.DiGraph:
    """水素結合ネットワークの有向グラフを生成する。

    無向グラフ（水素結合ネットワーク）から、各水素結合の方向（プロトンの向き）
    を決定して有向グラフを生成します。固定エッジで指定された方向は維持され、
    それ以外のエッジは双極子最適化アルゴリズムにより方向が決定されます。

    Args:
        graph: 拡大単位胞全体の無向グラフ
        depol_loop: 双極子最適化の反復回数
        lattice_sites: 格子サイト位置の配列
        fixedEdges: 固定エッジの有向グラフ

    Returns:
        nx.DiGraph: 水素結合の方向が決定された有向グラフ
    """
    for edge in fixedEdges.edges():
        _genice3_logger.debug(f"+ {edge=}")
    dg = genice_core.ice_graph(
        graph,
        vertexPositions=lattice_sites,
        isPeriodicBoundary=True,
        dipoleOptimizationCycles=depol_loop,
        fixedEdges=fixedEdges,
    )
    if not dg:
        raise ConfigurationError("Failed to generate a directed graph.")
    return dg


def orientations(
    lattice_sites: np.ndarray,
    digraph: nx.DiGraph,
    cell: np.ndarray,
    anions: Dict[int, str],
    cations: Dict[int, str],
) -> np.ndarray:
    """各水分子の配向行列を計算する。

    有向グラフで決定された水素結合の方向に基づいて、各水分子の
    配向（回転行列）を計算します。水分子のOHベクトルの方向から
    分子全体の配向を決定します。

    Args:
        lattice_sites: 格子サイト位置の配列
        digraph: 水素結合の方向が決定された有向グラフ
        cell: 拡大単位胞のセル行列
        anions: アニオン配置の辞書
        cations: カチオン配置の辞書

    Returns:
        np.ndarray: 各水分子の配向行列（Nx3x3配列、Nは水分子数）
    """
    return _assume_water_orientations(
        lattice_sites,
        digraph,
        cell,
        anions | cations,
    )


def cages(
    unitcell: UnitCell,
    replica_vectors: np.ndarray,
    replication_matrix: np.ndarray,
) -> CageSpecs:
    """ケージ位置とタイプを拡大単位胞全体に複製する。

    基本単位胞内で定義されたゲスト分子を配置するためのケージ（空隙）の
    位置とタイプを、拡大単位胞全体に複製します。

    Args:
        unitcell: 基本単位胞オブジェクト
        replica_vectors: レプリカベクトルの配列
        replication_matrix: 単位胞を複製するための3x3整数行列

    Returns:
        CageSpecs: 拡大単位胞全体でのケージ位置とタイプの情報
                   （unitcell.cagesがNoneの場合はNoneを返す）
    """
    if unitcell.cages is None:
        return None
    # ケージが存在しない場合（空の配列）の処理
    if len(unitcell.cages.positions) == 0:
        return CageSpecs(positions=np.array([]).reshape(0, 3), specs=[])
    repcagepos = replicate_positions(
        unitcell.cages.positions, replica_vectors, replication_matrix
    )
    num_cages_in_unitcell = len(unitcell.cages.positions)
    repcagespecs = [
        unitcell.cages.specs[i % num_cages_in_unitcell] for i in range(len(repcagepos))
    ]
    # unit cellのcagesと同じ構造。
    _genice3_logger.debug(f"{repcagepos=}, {repcagespecs=}")
    return CageSpecs(positions=repcagepos, specs=repcagespecs)


def cage_survey(
    cages: CageSpecs,
) -> str:
    """ケージ情報をJSON形式の文字列として返す。

    Args:
        cages: ケージ位置とタイプの情報

    Returns:
        str: ケージ情報をJSON形式でシリアライズした文字列
    """
    # JSONで返せばいい。ただ、graphなどは、もっと噛みくだいた形にしたいし、座標はnumpy arrayのままではJSONにならない。
    data = cages.to_json_capable_data()
    return json.dumps(data, indent=4)


# ============================================================================
# GenIce3クラス: DependencyEngineをラップ
# ============================================================================


class GenIce3:
    """GenIce3のメインクラス：リアクティブプロパティによる氷構造生成システム

    GenIce3は、依存関係エンジン（DependencyEngine）を使用して、必要な時に
    自動的に計算されるリアクティブプロパティを提供します。これにより、ユーザーは
    必要なプロパティにアクセスするだけで、そのプロパティに依存するすべての
    計算が自動的に実行されます。

    リアクティブプロパティの仕組み:
        - 各プロパティ（digraph, graph, lattice_sitesなど）は、アクセス時に
          必要に応じて自動的に計算されます
        - 依存関係は関数の引数名から自動的に推論されます
        - 一度計算されたプロパティはキャッシュされ、依存する入力が変更されるまで
          再利用されます
        - 入力プロパティ（unitcell, replication_matrix, depol_loopなど）が
          変更されると、それに依存するすべてのプロパティのキャッシュが自動的に
          クリアされます

    使用例:
        >>> genice = GenIce3(unitcell=my_unitcell)
        >>> digraph = genice.digraph  # 自動的に必要な計算が実行される
        >>> orientations = genice.orientations  # digraphに依存するため、digraphが先に計算される

    Attributes:
        digraph (nx.DiGraph): 水素結合の方向が決定された有向グラフ。
            無向グラフから双極子最適化アルゴリズムにより各水素結合の方向を決定した
            有向グラフです。固定エッジで指定された方向は維持され、それ以外のエッジは
            最適化により決定されます。このプロパティにアクセスすると、必要な依存関係
            （graph, lattice_sites, fixedEdgesなど）が自動的に計算されます。

        graph (nx.Graph): 水素結合ネットワークの無向グラフ。
            拡大単位胞全体の水分子間の水素結合ネットワークを表す無向グラフです。
            基本単位胞のグラフを拡大単位胞全体に複製して統合したものです。

        lattice_sites (np.ndarray): 格子サイト位置の配列。
            拡大単位胞内のすべての水分子の座標を表すNx3配列です（Nは水分子数）。
            基本単位胞内の水分子の座標を、単位胞の周期境界条件に従って拡大単位胞全体に
            複製したものです。

        orientations (np.ndarray): 各水分子の配向行列。
            各水分子の配向（回転行列）を表すNx3x3配列です（Nは水分子数）。
            有向グラフで決定された水素結合の方向に基づいて、各水分子のOHベクトルの
            方向から分子全体の配向を決定します。

        unitcell (UnitCell): 基本単位胞オブジェクト。
            氷構造の基本単位胞を表すオブジェクトです。格子構造、水分子の配置、
            水素結合ネットワーク、固定エッジなどの情報を含みます。

        replication_matrix (np.ndarray): 単位胞複製行列。
            基本単位胞をどのように積み重ねて拡大単位胞を作成するかを指定する3x3整数行列です。
            単位行列の場合は基本単位胞のみを使用します。

        depol_loop (int): 双極子最適化の反復回数。
            有向グラフを生成する際の双極子最適化アルゴリズムの反復回数です。
            値が大きいほどより最適化された構造が得られますが、計算時間も増加します。

        seed (int): 乱数シード。
            乱数生成器のシード値です。digraphの生成などで使用される乱数の初期化に使用されます。
            このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
            キャッシュが自動的にクリアされます。

        spot_anions (Dict[int, str]): 特定の格子サイト位置に配置するアニオンイオンの辞書。
            サイトインデックスからイオン名へのマッピングです。このプロパティを変更すると、
            それに依存するすべてのリアクティブプロパティのキャッシュが自動的にクリアされます。

        spot_cations (Dict[int, str]): 特定の格子サイト位置に配置するカチオンイオンの辞書。
            サイトインデックスからイオン名へのマッピングです。このプロパティを変更すると、
            それに依存するすべてのリアクティブプロパティのキャッシュが自動的にクリアされます。
    """

    # Class名でlog表示したい。
    logger = getLogger("GenIce3")

    # ユーザー向けAPIとして公開するpropertyのリスト
    # このリストに含まれるpropertyのみAPIドキュメントを作成する
    PUBLIC_API_PROPERTIES = [
        "digraph",
        "graph",
        "lattice_sites",
        "orientations",
        "unitcell",
        "replication_matrix",
        "depol_loop",
        "seed",
        "spot_anions",
        "spot_cations",
    ]

    def __init__(
        self,
        depol_loop: int = 1000,
        replication_matrix: np.ndarray = np.eye(3, dtype=int),
        seed: int = 1,
        spot_anions: Dict[int, str] = {},
        spot_cations: Dict[int, str] = {},
        **kwargs: Any,
    ) -> None:
        """GenIce3インスタンスを初期化する。

        Args:
            depol_loop: 双極子最適化の反復回数（デフォルト: 1000）
            replication_matrix: 単位胞複製行列（デフォルト: 単位行列）
            seed: 乱数シード（デフォルト: 1）
            spot_anions: 特定位置のアニオン配置（デフォルト: {}）
            spot_cations: 特定位置のカチオン配置（デフォルト: {}）
            **kwargs: その他のリアクティブプロパティ（unitcellなど）

        Raises:
            ConfigurationError: 無効なキーワード引数が指定された場合
        """
        # DependencyEngineインスタンスを作成
        self.engine = DependencyEngine()

        # Default値が必要なもの
        self.seed = seed  # reactive propertyとして設定（setterでnp.random.seed()も実行される）
        self.depol_loop = depol_loop
        self.replication_matrix = replication_matrix
        self.spot_anions = spot_anions
        self.spot_cations = spot_cations

        # タスクを登録（モジュールレベルの関数を登録）
        self._register_tasks()

        # Default値が不要なもの
        for key in self.list_settable_reactive_properties():
            if key in kwargs:
                setattr(self, key, kwargs.pop(key))
        if kwargs:
            raise ConfigurationError(f"Invalid keyword arguments: {kwargs}.")

    def _register_tasks(self):
        """DependencyEngineにモジュールレベルのタスク関数を登録する"""
        engine = self.engine
        # モジュールから関数を取得して登録
        import sys

        current_module = sys.modules[__name__]

        task_names = [
            "cell",
            "replica_vectors",
            "replica_vector_labels",
            "graph",
            "lattice_sites",
            "anions",
            "cations",
            "site_occupants",
            "fixedEdges",
            "digraph",
            "orientations",
            "cages",
            "cage_survey",
        ]

        for name in task_names:
            func = getattr(current_module, name)
            engine.task(func)

    # spot_anions
    @property
    def spot_anions(self):
        """特定の格子サイト位置に配置するアニオンイオンの辞書。

        Returns:
            Dict[int, str]: サイトインデックスからイオン名へのマッピング
        """
        if not hasattr(self, "_spot_anions"):
            self._spot_anions = {}
        return self._spot_anions

    @spot_anions.setter
    def spot_anions(self, spot_anions):
        """特定の格子サイト位置に配置するアニオンイオンを設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            spot_anions: サイトインデックスからイオン名へのマッピング
        """
        self._spot_anions = spot_anions
        self.logger.debug(f"  {spot_anions=}")
        self.engine.cache.clear()

    # spot_cations
    @property
    def spot_cations(self):
        """特定の格子サイト位置に配置するカチオンイオンの辞書。

        Returns:
            Dict[int, str]: サイトインデックスからイオン名へのマッピング
        """
        if not hasattr(self, "_spot_cations"):
            self._spot_cations = {}
        return self._spot_cations

    @spot_cations.setter
    def spot_cations(self, spot_cations):
        """特定の格子サイト位置に配置するカチオンイオンを設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            spot_cations: サイトインデックスからイオン名へのマッピング
        """
        self._spot_cations = spot_cations
        self.logger.debug(f"  {spot_cations=}")
        self.engine.cache.clear()

    @property
    def seed(self):
        """乱数シード。

        Returns:
            int: 乱数シード値
        """
        return self._seed

    @seed.setter
    def seed(self, seed):
        """乱数シードを設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            seed: 乱数シード値
        """
        self._seed = seed
        np.random.seed(seed)
        self.logger.debug(f"  {seed=}")
        # キャッシュをクリア（seedに依存するすべてのタスクを再計算させる）
        self.engine.cache.clear()

    @property
    def depol_loop(self):
        """双極子最適化の反復回数。

        Returns:
            int: 双極子最適化アルゴリズムの反復回数
        """
        return self._depol_loop

    @depol_loop.setter
    def depol_loop(self, depol_loop):
        """双極子最適化の反復回数を設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            depol_loop: 双極子最適化アルゴリズムの反復回数
        """
        self._depol_loop = depol_loop
        self.logger.debug(f"  {depol_loop=}")
        # キャッシュをクリア（depol_loopに依存するすべてのタスクを再計算させる）
        self.engine.cache.clear()

    @property
    def unitcell(self):
        """基本単位胞オブジェクト。

        Returns:
            UnitCell: 基本単位胞オブジェクト

        Raises:
            ConfigurationError: 単位胞が設定されていない場合
        """
        if not hasattr(self, "_unitcell") or self._unitcell is None:
            raise ConfigurationError("Unitcell is not set.")
        return self._unitcell

    @unitcell.setter
    def unitcell(self, unitcell):
        """基本単位胞オブジェクトを設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            unitcell: 基本単位胞オブジェクト
        """
        self._unitcell = unitcell
        self.logger.debug(f"  {unitcell=}")
        self.logger.debug(f"  {unitcell.waters=}")
        self.logger.debug(f"  {unitcell.graph=}")
        self.logger.debug(f"  {unitcell.fixed=}")

        a, b, c, A, B, C = cellshape(self.replication_matrix @ self.unitcell.cell)
        self.logger.debug("  Reshaped cell:")
        self.logger.debug(f"    {a=:.4f}, {b=:.4f}, {c=:.4f}")
        self.logger.debug(f"    {A=:.3f}, {B=:.3f}, {C=:.3f}")
        # キャッシュをクリア（unitcellに依存するすべてのタスクを再計算させる）
        self.engine.cache.clear()

    @property
    def replication_matrix(self):
        """単位胞を複製するための3x3整数行列。

        この行列により、基本単位胞をどのように積み重ねて拡大単位胞を
        作成するかを指定します。

        Returns:
            np.ndarray: 3x3整数行列
        """
        return self._replication_matrix

    @replication_matrix.setter
    def replication_matrix(self, replication_matrix):
        """単位胞複製行列を設定する。

        このプロパティを変更すると、それに依存するすべてのリアクティブプロパティの
        キャッシュが自動的にクリアされます。

        Args:
            replication_matrix: 3x3整数行列
        """
        self._replication_matrix = replication_matrix
        i, j, k = np.array(replication_matrix)
        self.logger.debug(f"    {i=}")
        self.logger.debug(f"    {j=}")
        self.logger.debug(f"    {k=}")
        # キャッシュをクリア（replication_matrixに依存するすべてのタスクを再計算させる）
        self.engine.cache.clear()

    def _get_inputs(self) -> Dict[str, Any]:
        """engine.resolve()に渡すinputs辞書を取得"""
        return {
            "unitcell": self.unitcell,
            "replication_matrix": self.replication_matrix,
            "depol_loop": self.depol_loop,
            "spot_anions": self.spot_anions,
            "spot_cations": self.spot_cations,
        }

    def __getattr__(self, name: str):
        """リアクティブプロパティへのアクセスを自動解決する。

        このメソッドにより、DependencyEngineに登録されたタスク関数が
        プロパティとしてアクセス可能になります。プロパティにアクセスすると、
        依存関係が自動的に解決され、必要な計算が実行されます。

        例:
            >>> genice = GenIce3(unitcell=my_unitcell)
            >>> digraph = genice.digraph  # この時点でdigraphとその依存関係が計算される

        Args:
            name: アクセスするプロパティ名（タスク関数名）

        Returns:
            計算されたプロパティの値

        Raises:
            AttributeError: 指定された名前のプロパティが存在しない場合
            ConfigurationError: 必要な入力が設定されていない場合
        """
        # DependencyEngineに登録されているタスクかチェック
        if name in self.engine.registry:
            return self.engine.resolve(name, self._get_inputs())

        # それ以外は通常のAttributeErrorを発生
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    def water_molecules(self, water_model: Molecule) -> Dict[int, Molecule]:
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

    def guest_molecules(
        self, guests: Dict[str, List[GuestSpec]], spot_guests: Dict[int, Molecule]
    ) -> List[Molecule]:

        # self.cagesがNoneの場合は、そもそもケージの調査をしていない可能性がある。
        # ケージの調査は本来はunitcellで行うべきことだが、そこまで依存関係で辿れない。

        # if self.cages is None:
        #     #
        #     if len(guests) > 0:
        #         raise ConfigurationError("Cages are not defined.")
        #     if len(spot_guests) > 0:
        #         raise ConfigurationError("Spot guests are not defined.")
        #     return []

        all_labels = [spec.label for spec in self.cages.specs]
        # guest_specで指定されているのに存在しない種類のケージはエラーにする。
        for label in guests:
            if label not in all_labels:
                raise ConfigurationError(f"Cage type {label} is not defined.")

        randoms = np.random.random(len(self.cages.positions))

        mols = []
        for pos, spec, probability in zip(
            self.cages.positions, self.cages.specs, randoms
        ):
            accum = 0.0
            if spec.label in guests:
                for guest_spec in guests[spec.label]:
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

        # spot guestの配置
        for cage_index, guest in spot_guests.items():
            mols.append(
                Molecule(
                    name=guest.name,
                    sites=guest.sites + self.cages.positions[cage_index] @ self.cell,
                    labels=guest.labels,
                    is_water=False,
                )
            )
        return mols

    def substitutional_ions(self) -> Dict[int, Molecule]:
        # 将来は分子イオン(H3O+など)を置く可能性もあることに留意。
        ions: Dict[int, Molecule] = {}
        # ならべかえはここではしない。formatterにまかせる。
        for label, molecule in self.anions.items():
            ions[label] = Molecule(
                name=molecule,
                sites=[self.lattice_sites[label] @ self.cell],
                labels=[molecule],
                is_water=False,
            )
        for label, molecule in self.cations.items():
            ions[label] = Molecule(
                name=molecule,
                sites=[self.lattice_sites[label] @ self.cell],
                labels=[molecule],
                is_water=False,
            )
        return ions

    def dope_anions(self, anions: Dict[int, Molecule]):
        for label, molecule in anions.items():
            self.anions[label] = molecule

    def dope_cations(self, cations: Dict[int, Molecule]):
        for label, molecule in cations.items():
            self.cations[label] = molecule

    # def get_atomic_structure(
    #     self,
    #     water_model: Optional[Molecule] = None,
    #     guests: Optional[Dict[str, List[GuestSpec]]] = None,
    #     spot_guests: Optional[Dict[int, Molecule]] = None,
    # ) -> AtomicStructure:
    #     """
    #     原子構造データを統合的に取得する。
    #     exporterプラグインが使用するための統一インターフェース。

    #     Args:
    #         water_model: 水分子モデル（デフォルト: FourSiteWater()）
    #         guests: ケージタイプごとのゲスト分子の指定（デフォルト: {}）
    #         spot_guests: 特定ケージ位置へのゲスト分子の指定（デフォルト: {}）

    #     Returns:
    #         AtomicStructure: 水分子、ゲスト分子、イオン、セル行列を含む統合データ構造
    #     """

    #     return AtomicStructure(
    #         waters=self.water_molecules(water_model=water_model),
    #         guests=self.guest_molecules(
    #             guests=guests or {}, spot_guests=spot_guests or {}
    #         ),
    #         ions=self.substitutional_ions(),
    #         cell=self.cell,
    #     )

    @classmethod
    def get_public_api_properties(cls):
        """
        ユーザー向けAPIとして公開されているpropertyのリストを返す。

        Returns:
            list: 公開API property名のリスト
        """
        return cls.PUBLIC_API_PROPERTIES.copy()

    def list_all_reactive_properties(self):
        """
        すべてのリアクティブプロパティ（DependencyEngineに登録されたタスク）を列挙する。

        Returns:
            dict: property名をキー、タスク関数を値とする辞書
        """
        return self.engine.registry.copy()

    def list_public_reactive_properties(self):
        """
        ユーザー向けAPIとして公開されているリアクティブプロパティのみを列挙する。

        Returns:
            dict: property名をキー、タスク関数を値とする辞書
        """
        all_reactive = self.list_all_reactive_properties()
        public_names = set(self.PUBLIC_API_PROPERTIES)
        return {
            name: func for name, func in all_reactive.items() if name in public_names
        }

    @classmethod
    def list_settable_reactive_properties(cls):
        """setterを持つリアクティブな変数を列挙する"""
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
