from dataclasses import dataclass
from typing import Dict, Set, List, Tuple
from logging import getLogger
from collections import defaultdict

import numpy as np
import networkx as nx

from genice2.decorators import banner, timeit
from genice2 import ConfigurationError


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


@dataclass
class Stage2Output:
    """Stage2の出力データ"""

    graph: nx.Graph  # ネットワークトポロジー
    dopants: Dict[int, str]  # ドーパント情報
    fixed_edges: nx.DiGraph  # 固定エッジ
    groups: Dict[int, Dict[int, str]]  # グループ情報
    filled_cages: Set[int]  # 埋められたケージ


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
        self.anions1 = anions
        self.cations1 = cations
        # self.dopants1 = anions | cations
        self.groups1 = groups1
        self.cagepos1 = cagepos1

    @timeit
    @banner
    def execute(self) -> Stage2Output:
        """Makes a random graph and replicates it."""
        logger = getLogger()

        # # Some edges are directed when ions are doped.
        # if self.doping_hook_function is not None:
        #     self.doping_hook_function(self)  # may be defined in the plugin

        # Replicate the dopants in the unit cell
        # dopants = replicate_labeldict(
        #     self.dopants1, len(self.waters1), self.replica_vectors
        # )
        anions = replicate_labeldict(
            self.anions1, len(self.waters1), self.replica_vectors
        )
        cations = replicate_labeldict(
            self.cations1, len(self.waters1), self.replica_vectors
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

        # Dope ions by options.
        # replicateされた構造にドープする。
        dopants = anions | cations
        if len(anions) > 0:
            logger.info(f"  Anionize: {anions}.")
            for site, _ in anions.items():
                for nei in graph[site]:
                    if fixed_edges.has_edge(site, nei):
                        raise ConfigurationError(
                            f"Impossible to dope an anion at {site}."
                        )
                    fixed_edges.add_edge(nei, site)
        if len(cations) > 0:
            logger.info(f"  Cationize: {cations}.")
            for site, _ in cations.items():
                for nei in graph[site]:
                    if fixed_edges.has_edge(nei, site):
                        raise ConfigurationError(
                            f"Impossible to dope a cation at {site}."
                        )
                    fixed_edges.add_edge(site, nei)

        # Count bonds
        nfixed = fixed_edges.number_of_edges()
        num_hb_disorder = graph.number_of_edges() - nfixed
        logger.info(f"  Number of pre-oriented hydrogen bonds: {nfixed}")
        logger.info(f"  Number of unoriented hydrogen bonds: {num_hb_disorder}")
        logger.info(f"  Total number of hydrogen bonds: {nfixed + num_hb_disorder}")

        # return Stage2Output(dopants, groups, filled_cages, graph, fixedEdges)
        return Stage2Output(
            graph=graph,
            dopants=dopants,
            fixed_edges=fixed_edges,
            groups=groups,
            filled_cages=filled_cages,
        )
