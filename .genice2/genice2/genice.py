"""
GenIce class
"""

import itertools as it
from collections import defaultdict
from logging import getLogger
from typing import Type, Union

import numpy as np
import pairlist as pl
import networkx as nx

from genice2 import cage

# from genice2 import digraph_unused as dg
from genice2.cell import Cell, rel_wrap, cellshape
from genice2.decorators import banner, timeit
from genice2.formats import Format
from genice2.lattices import Lattice
import genice_core

# A virtual monatomic molecule
from genice2.molecules import Molecule, arrange, monatom, one
from genice2.plugin import Group, safe_import
from genice2.valueparser import (
    parse_cages,
    parse_pairs,
    plugin_option_parser,
    put_in_array,
)


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


def orientations(coord, graph, cell, immutables: set):
    """
    Does not work when two OHs are colinear
    """

    logger = getLogger()
    # just for a test of pure water
    assert len(coord) == graph.number_of_nodes(), (len(coord), graph.number_of_nodes())
    logger.info(f"{immutables} immutables")
    # 通常の氷であればアルゴリズムを高速化できる。

    nnode = len(list(graph))
    neis = np.zeros([nnode, 2], dtype=int)

    # 仮想ノード用の配列。第0要素は実際には第nnode要素を表す。
    extended_coord = []

    # v0 = np.zeros([nnode, 3])
    # v1 = np.zeros([nnode, 3])
    for node in graph:
        if node in immutables:
            h1 = np.array([0.0, 1, 1]) / (2**0.5)
            h2 = np.array([0.0, -1, 1]) / (2**0.5)
            r1 = cell.abs2rel(h1)
            r2 = cell.abs2rel(h2)
            # 仮想ノードにさしかえる
            neis[node] = [nnode + len(extended_coord), nnode + len(extended_coord) + 1]
            extended_coord += [coord[node] + r1, coord[node] + r2]
            continue
        succ = list(graph.successors(node))
        if len(succ) < 2:
            vsucc = cell.rel2abs(rel_wrap(coord[succ] - coord[node]))
            pred = list(graph.predecessors(node))
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


@timeit
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


class GenIce:
    """
    The core of GenIce.

    lat:     An instance of a Lattice class.
    density: Density of target ice in g/cm3.
    rep:     Repetition of the cell. A tuple of three integers.
    anions, cations:
             The locations of monovalent anions and cations that replace
             the water molecules.
    spot_guests:
             The locations of guest molecules that occupy the specified cages.
    spot_group: (EXPERIMENTAL)
             The locations of functional groups that occupy the cages.
    as_is:   Avoids shuffling of the orientations of water molecules.
    signature: A text that is inserted in the output.
    reshape: reshape matrix.
    """

    @timeit
    @banner
    def __init__(
        self,
        lat: Type[Lattice],
        signature: str = "",
        density: float = 0,
        rep=None,
        reshape=np.eye(3, dtype=int),
        cations: dict = {},
        anions: dict = {},
        spot_guests: dict = {},
        spot_groups: dict = {},
        asis: bool = False,
        shift=(0.0, 0.0, 0.0),
        # seed=1000,
    ):
        """
        Constructor of GenIce.

        Arguments:
            lat:        The ice lattice.
            signature:  A string for a signature.
            density:    Target density.
            rep:        Cell repetitions.
            cations:    Labels of water molecules that are replaced by the cations.
            anions:
            spot_guests:Labels of cages in which a guest is placed.
            spot_groups:Labels of cages in which a group is placed.
            asis:       Do not modify the orientations of the hydrogen bonds.
            shift:      A fractional value to be added to the positions.
            reshape: reshape matrix.
        """

        # このconstructorが大きすぎ、複雑すぎ。
        # self変数も多すぎ。
        logger = getLogger()
        self.asis = asis
        self.cations = cations
        self.anions = anions
        self.spot_guests = spot_guests
        self.spot_groups = spot_groups

        # 変数名に1がついているのはunit cell

        # ================================================================
        # cell: cell dimension
        #   see parse_cell for syntax.
        #
        self.cell1 = Cell(lat.cell)

        # Reshaping matrix (Must be integers)
        # for now they are hardcoded.  It will be given as options for the plugin in the future.
        # ijk = np.array([[1, 1, 1], [1, -1, 0], [1, 1, -2]])
        # ijk = np.array([[2, 0, 1], [0, 1, 0], [-1, 0, 2]])
        # ijk = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 2]])
        # ijk = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])

        # CLIからはrepは与えられない。API経由で設定される可能性はある。
        if rep is not None:
            logger.warning("rep for GenIce() is deprecated. Use reshape instead.")
            # 直方体レプリカの順序指定。
            self.replica_vectors = np.array(
                [
                    (x, y, z)
                    for x in range(rep[0])
                    for y in range(rep[1])
                    for z in range(rep[2])
                ]
            )
            # セルの複製後の大セルの形
            self.reshape_matrix = np.diag(rep)
        else:
            logger.info("  Reshaping the unit cell.")

            self.reshape_matrix = reshape

            i, j, k = np.array(reshape)
            logger.info(f"    i:{i}")
            logger.info(f"    j:{j}")
            logger.info(f"    k:{k}")

            a, b, c, A, B, C = cellshape(reshape @ self.cell1.mat)
            logger.info("  Reshaped cell:")
            logger.info(f"    a,b,c = {a}, {b}, {c}")
            logger.info(f"    A,B,C = {A}, {B}, {C}")
            abc = reshape @ self.cell1.mat

            # 単位胞を並べた格子を、大セルで切りとった時に、どの範囲の単位胞がかするか
            corners = np.array(
                [a * i + b * j + c * k for a in (0, 1) for b in (0, 1) for c in (0, 1)]
            )

            mins = np.min(corners, axis=0)
            maxs = np.max(corners, axis=0)

            # 整数行列の逆行列に行列式をかけたものは整数行列になる。
            det = np.linalg.det(reshape)
            if det < 0:
                det = -det
            det = np.floor(det + 0.5).astype(int)
            invdet = np.floor(np.linalg.inv(reshape) * det + 0.5).astype(int)

            vecs = set()
            # かする単位胞のうち
            for a in range(mins[0], maxs[0] + 1):
                logger.debug(a)
                for b in range(mins[1], maxs[1] + 1):
                    for c in range(mins[2], maxs[2] + 1):
                        # 単位胞の位置を、
                        abc = np.array([a, b, c])
                        # 大セルにおさめ
                        rep = grandcell_wrap(abc, reshape, invdet, det).astype(int)
                        # 記録する
                        if tuple(rep) not in vecs:
                            vecs.add(tuple(rep))

            self.replica_vectors = np.array(list(vecs))

            # 大セルの大きさは、単位胞の整数倍でなければいけない。
            vol = abs(np.linalg.det(reshape))
            assert np.allclose(vol, len(vecs)), (vol, vecs)

            # 大セルの形
            self.reshape_matrix = reshape

        # レプリカ単位胞には0から順に番号ラベルがついている。replica_vector_labelsは単位胞の位置をラベルに変換する
        self.replica_vector_labels = {
            tuple(xyz): i for i, xyz in enumerate(self.replica_vectors)
        }

        # Show the document of the module
        try:
            self.doc = lat.__doc__.splitlines()
        except BaseException:
            self.doc = []

        if len(signature) > 0:
            self.doc.append("")
            self.doc.append(signature)

        for line in self.doc:
            logger.info("  " + line)

        # ================================================================
        # waters: positions of water molecules
        #
        self.waters1 = put_in_array(lat.waters)
        logger.debug(f"Waters: {len(self.waters1)}")
        self.waters1 = self.waters1.reshape((-1, 3))

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative", i.e. in fractional coordinates
        #
        if lat.coord == "absolute":
            self.waters1 = self.cell1.abs2rel(self.waters1)

        # shift of the origin
        self.waters1 = np.array(self.waters1) + np.array(shift)
        # fractional coordinate between [0, 1)
        self.waters1 -= np.floor(self.waters1)

        # ================================================================
        # pairs: specify the pairs of molecules that are connected.
        #   Bond orientation will be shuffled later
        #   unless it is "fixed".
        #
        self.pairs1 = None

        try:
            self.pairs1 = parse_pairs(lat.pairs)
        except AttributeError:
            logger.info("HB connectivity is not defined.")

        # ================================================================
        # bondlen: specify the bond length threshold.
        #   This is used when "pairs" are not specified.
        #   It is applied to the original positions of molecules (before density setting).
        #
        nmol = self.waters1.shape[0]  # nmol in a unit cell
        volume = self.cell1.volume()  # volume of a unit cell in nm**3
        self.bondlen = None

        try:
            self.bondlen = lat.bondlen
            logger.info(f"Bond length (specified): {self.bondlen}")
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            # very rough estimate of the pair distances assuming that the particles distribute homogeneously.
            rc = (volume / nmol) ** (1 / 3) * 1.5
            p = pl.pairs_iter(
                self.waters1, maxdist=rc, cell=self.cell1.mat, distance=False
            )
            self.bondlen = 1.1 * shortest_distance(self.waters1, self.cell1, pairs=p)
            logger.info(f"Bond length (estim.): {self.bondlen}")

        # Set density
        mass = 18  # water
        NB = 6.022e23
        density0 = mass * nmol / (NB * volume * 1e-21)

        if density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                logger.info(
                    "Density is not specified. Assume the density from lattice."
                )
                dmin = shortest_distance(self.waters1, self.cell1)
                logger.info(
                    f"Closest pair distance: {dmin} (should be around 0.276 nm)"
                )
                self.density = density0 / (0.276 / dmin) ** 3
                # self.density = density0
        else:
            self.density = density

        logger.info(f"Target Density: {self.density}")
        logger.info(f"Original Density: {density0}")

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density) ** (1.0 / 3.0)
        self.cell1.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        logger.info(f"Bond length (scaled, nm): {self.bondlen}")

        # ================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        self.cagepos1 = None
        self.cagetype1 = None

        if "cages" in lat.__dict__:
            self.cagepos1, self.cagetype1 = parse_cages(lat.cages)
            logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            # pre-parsed data
            self.cagepos1, self.cagetype1 = np.array(lat.cagepos), lat.cagetype
        if self.cagepos1 is not None:
            self.cagepos1 = np.array(self.cagepos1) + np.array(shift)
            self.cagepos1 -= np.floor(self.cagepos1)

        # ================================================================
        # fixed: specify the bonds whose directions are fixed.
        #   you can specify them in pairs at a time.
        #   You can also leave it undefined.
        #
        self.fixed1 = []
        try:
            self.fixed1 = parse_pairs(lat.fixed)
            logger.info("Orientations of some edges are fixed.")
        except AttributeError:
            pass

        if "dopeIonsToUnitCell" in lat.__dict__:
            self.dopeIonsToUnitCell = lat.dopeIonsToUnitCell
        else:
            self.dopeIonsToUnitCell = None
        self.dopants1 = set()

        # self.immutables1 = set(self.dopants1)
        # logger.info(f"{self.dopants1} dopants1")

        # if asis, make pairs to be fixed.
        if self.asis and len(self.fixed1) == 0:
            self.fixed1 = self.pairs1

        # filled cages
        self.filled_cages = set()

        # groups info
        self.groups1 = defaultdict(dict)

    def generate_ice(
        self,
        formatter: Type[Format],
        water: Union[Type[Molecule], None] = None,
        guests={},
        depol="strict",
        noise=0.0,
        assess_cages=False,
    ):
        """
        Generate an ice structure and dump it with the aid of a formatter plugin.

            formatter: genice2.format.Format() class
            water:     genice2.molecules.Molecule() class
            assess_cages:   Cages will be assessed on the fly instead of
                        pre-specified in the lattice plugin.
        """

        logger = getLogger()

        # in old syntax, the arguments water and formatter were mandatory, but
        # in new syntax, water is optional and their order is exchanged.
        # therefore i prepare a backward compatibility.
        from genice2.molecules import Molecule

        if isinstance(formatter, Molecule):
            formatter, water = water, formatter
            logger.warn(
                "generate_ice(water, formatter) is deprecated. "
                "New syntax is: generate_ice(formatter, water=water)."
            )

        def Stages():
            hooks = formatter.hooks()
            maxstage = max(0, *hooks.keys())

            if 0 in hooks:
                abort = hooks[0](self)
                if maxstage < 1 or abort:
                    return

            self.Stage1(noise, assess_cages=assess_cages)

            if 1 in hooks:
                abort = hooks[1](self)
                if maxstage < 2 or abort:
                    return

            self.Stage2()

            # Count bonds
            nfixed = self.fixedEdges.number_of_edges()
            num_hb_disorder = self.graph.number_of_edges() - nfixed
            logger.info(f"  Number of pre-oriented hydrogen bonds: {nfixed}")
            logger.info(f"  Number of unoriented hydrogen bonds: {num_hb_disorder}")
            logger.info(
                "  Number of hydrogen bonds: {0} (regular num: {1})".format(
                    nfixed + num_hb_disorder, len(self.reppositions) * 2
                )
            )

            if 2 in hooks:
                abort = hooks[2](self)
                if maxstage < 3 or abort:
                    return

            # GenIce-core
            self.Stage34E(depol=depol)

            if 3 in hooks:
                abort = hooks[3](self)
                if maxstage < 4 or abort:
                    return

            if 4 in hooks:
                abort = hooks[4](self)
                if maxstage < 5 or abort:
                    return

            self.Stage5()

            if 5 in hooks:
                abort = hooks[5](self)
                if maxstage < 6 or abort:
                    return

            self.Stage6(water)

            if 6 in hooks:
                abort = hooks[6](self)
                if maxstage < 7 or abort:
                    return

            self.Stage7(guests)

            if 7 in hooks:
                hooks[7](self)

        abort = Stages()
        if not abort:
            return formatter.dump()

    @timeit
    @banner
    def Stage1(self, noise=0.0, assess_cages=False):
        """
        Replicate water molecules to make a repeated cell.

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        repcagetype: replicated cage types array
        repcagepos:  replicated cage positions (CoM, relative)
        cagetypes:   set of cage types
        """

        logger = getLogger()
        self.reppositions = replicate_positions(
            self.waters1, self.replica_vectors, self.reshape_matrix
        )

        # This must be done before the replication of the cell.
        logger.info("  Number of water molecules: {0}".format(len(self.reppositions)))

        if self.pairs1 is None:
            self.pairs1 = self.prepare_pairs()

        self.graph1 = nx.Graph(self.pairs1)
        logger.info(f"  Number of water nodes: {self.graph1.number_of_nodes()}")

        # scale the cell
        self.repcell = Cell(self.reshape_matrix @ self.cell1.mat)

        if noise > 0.0:
            logger.info(f"  Add noise: {noise}.")
            perturb = np.random.normal(
                loc=0.0,
                scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                size=self.reppositions.shape,
            )
            self.reppositions += self.repcell.abs2rel(perturb)

        if assess_cages:
            logger.info("  Assessing the cages...")
            self.cagepos1, self.cagetype1 = cage.assess_cages(self.graph1, self.waters1)
            logger.info("  Done assessment.")

        if self.cagepos1 is not None and self.cagepos1.shape[0] > 0:
            logger.info("  Hints:")
            self.repcagepos = replicate_positions(
                self.cagepos1, self.replica_vectors, self.reshape_matrix
            )
            nrepcages = self.repcagepos.shape[0]
            self.repcagetype = [
                self.cagetype1[i % len(self.cagetype1)] for i in range(nrepcages)
            ]
            self.cagetypes = defaultdict(set)

            for i, typ in enumerate(self.repcagetype):
                self.cagetypes[typ].add(i)

            # INFO for cage types
            logger.info("    Cage types: {0}".format(list(self.cagetypes)))

            for typ, cages in self.cagetypes.items():
                logger.info(f"    Cage type {typ}: {cages}")
            # Up here move to stage 1.
        else:
            self.repcagetype = None
            self.repcagepos = None
            self.cagetypes = None

    @timeit
    @banner
    def Stage2(self):
        """
        Make a random graph and replicate.

        Provided variables:
        dopants:
        groups:  replicated positions of the chemical groups (CoM, relative)
        filled_cages:
        graph:   replicated network topology (bond orientation may be random)
        """

        logger = getLogger()

        # Some edges are directed when ions are doped.
        if self.dopeIonsToUnitCell is not None:
            self.dopeIonsToUnitCell(self)  # may be defined in the plugin

        # Replicate the dopants in the unit cell
        self.dopants = replicate_labeldict(
            self.dopants1, len(self.waters1), self.replica_vectors
        )

        self.groups = replicate_groups(
            self.groups1, self.waters1, self.cagepos1, self.replica_vectors
        )

        # self.groups_info(self.groups)
        for _, cages in self.groups.items():
            self.filled_cages |= set(cages)

        # logger.info(("filled",self.filled_cages))
        # Replicate the graph
        self.graph, self.fixedEdges = replicate_graph(
            self.graph1,
            self.waters1,
            self.replica_vectors,
            self.fixed1,
            self.replica_vector_labels,
            self.reshape_matrix,
        )
        # replicate "immutables" == dopants list in the graph
        # self.repimmutables = replicate_labels(
        #     self.immutables1, self.waters1.shape[0], self.rep
        # )

        self.repimmutables = set()
        # Dope ions by options.
        if len(self.anions) > 0:
            logger.info(f"  Anionize: {self.anions}.")

            for site, name in self.anions.items():
                # self.graph.anionize(site)
                for nei in self.graph[site]:
                    self.fixedEdges.add_edge(nei, site)
                self.dopants[site] = name
                self.repimmutables.add(site)

        if len(self.cations) > 0:
            logger.info(f"  Cationize: {self.cations}.")

            for site, name in self.cations.items():
                # self.graph.cationize(site)
                for nei in self.graph[site]:
                    self.fixedEdges.add_edge(site, nei)
                self.dopants[site] = name
                self.repimmutables.add(site)

    @timeit
    @banner
    def Stage34E(self, depol: str = "none"):
        """
        Make a directed ice graph from an undirected ice graph.
        """

        logger = getLogger()

        if depol == "none":
            iter = 0
        else:
            iter = 1000
        self.digraph = genice_core.ice_graph(
            self.graph.to_undirected(),
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=iter,
            fixedEdges=self.fixedEdges,
        )

        # self.digraph = nx.compose(d, self.fixedEdges)

        logger.debug("  Depolarized in Stage34E.")

    @timeit
    @banner
    def Stage5(self):
        """
        Prepare orientations for rigid molecules.

        Provided variables:
        reppositions: molecular positions.
        rotmatrices:  rotation matrices for water molecules
        """

        logger = getLogger()

        # determine the orientations of the water molecules based on edge
        # directions.
        self.rotmatrices = orientations(
            self.reppositions, self.digraph, self.repcell, self.repimmutables
        )

        # Activate it.
        # logger.info("The network is not specified.  Water molecules will be orinented randomly.")
        # rotmatrices = [rigid.rand_rotation_matrix() for pos in positions]

    @timeit
    @banner
    def Stage6(self, water):
        """
        Arrange atoms of water and replacements

        Provided variables:
        atoms: atomic positions of water molecules. (absolute)
        """

        logger = getLogger()

        # assert audit_name(water_type), "Dubious water name: {0}".format(water_type)
        # water = importlib.import_module("genice.molecules."+water_type)
        # water = safe_import("molecule", water_type)

        try:
            mdoc = water.__doc__.splitlines()
        except BaseException:
            mdoc = []

        for line in mdoc:
            logger.info("  " + line)

        self.universe = []
        self.universe.append(
            arrange(
                self.reppositions,
                self.repcell,
                self.rotmatrices,
                water,
                immutables=set(self.dopants),
            )
        )

    @timeit
    @banner
    def Stage7(self, guests):
        """
        Arrange guest atoms.

        Provided variables:
        atoms: atomic positions of all molecules.
        """

        logger = getLogger()

        if self.cagepos1 is not None:
            # the cages around the dopants.
            dopants_neighbors = dopants_info(
                self.dopants, self.reppositions, self.repcagepos, self.repcell
            )

            # put the (one-off) groups
            if len(self.spot_groups) > 0:
                # process the -H option
                for cage, group_to in self.spot_groups.items():
                    group, root = group_to.split(":")
                    self.add_group(cage, group, int(root))

            molecules = defaultdict(list)

            if len(self.spot_guests) > 0:
                # process the -G option
                for cage, molec in self.spot_guests.items():
                    molecules[molec].append(cage)
                    self.filled_cages.add(cage)

            # process the -g option
            for cagetype, contents in guests.items():
                if cagetype not in self.cagetypes:
                    logger.info(f"Nonexistent cage type: {cagetype}")
                    continue
                resident = dict()
                rooms = list(self.cagetypes[cagetype] - self.filled_cages)

                for room in rooms:
                    resident[room] = None

                vacant = len(rooms)

                for molec, frac in contents.items():
                    nmolec = int(frac * len(rooms) + 0.5)
                    vacant -= nmolec
                    assert vacant >= 0, "Too many guests."
                    remain = nmolec
                    movedin = []

                    while remain > 0:
                        r = np.random.randint(0, len(rooms))
                        room = rooms[r]

                        if resident[room] is None:
                            resident[room] = molec
                            molecules[molec].append(room)
                            movedin.append(room)
                            remain -= 1

            # Now ge got the address book of the molecules.
            if len(molecules):
                logger.info("  Summary of guest placements:")
                guests_info(self.cagetypes, molecules)

            if len(self.spot_groups) > 0:
                logger.info("  Summary of groups:")
                groups_info(self.groups)

            # semi-guests
            for root, cages in self.groups.items():
                assert root in self.dopants
                name = self.dopants[root]
                molname = f"G{root}"
                pos = self.reppositions[root]
                rot = self.rotmatrices[root]
                self.universe.append(monatom(pos, self.repcell, name))
                del self.dopants[root]  # processed.
                logger.debug((root, cages, name, molname, pos, rot))

                for cage, group in cages.items():
                    # assert group in self.groups_placer
                    assert cage in dopants_neighbors[root]
                    cage_center = self.repcagepos[cage]
                    self.universe.append(
                        Group(group).arrange_atoms(
                            cage_center, pos, self.repcell, molname, origin_atom=name
                        )
                    )

            # molecular guests
            for molec, cages in molecules.items():
                guest_type, guest_options = plugin_option_parser(molec)
                logger.debug(f"Guest type: {guest_type}")
                gmol = safe_import("molecule", guest_type).Molecule(**guest_options)

                try:
                    mdoc = gmol.__doc__.splitlines()
                except BaseException:
                    mdoc = []
                for line in mdoc:
                    logger.info("  " + line)
                cage_center = [self.repcagepos[i] for i in cages]
                cmat = [np.identity(3) for i in cages]
                self.universe.append(arrange(cage_center, self.repcell, cmat, gmol))

        # Assume the dopant is monatomic and replaces one water molecule
        atomset = defaultdict(set)
        for label, name in self.dopants.items():
            atomset[name].add(label)

        for name, labels in atomset.items():
            pos = [self.reppositions[i] for i in sorted(labels)]
            rot = [self.rotmatrices[i] for i in sorted(labels)]
            oneatom = one.Molecule(label=name)
            self.universe.append(arrange(pos, self.repcell, rot, oneatom))

    def prepare_pairs(self):
        logger = getLogger()

        logger.info("  Pairs are not given explicitly.")
        logger.info("  Estimating the bonds according to the pair distances.")

        logger.debug(f"Bondlen: {self.bondlen}")
        # make bonded pairs according to the pair distance.
        # make before replicating them.
        return [
            (i, j)
            for i, j in pl.pairs_iter(
                self.waters1, self.bondlen, self.cell1.mat, distance=False
            )
        ]

    def add_group(self, cage, group, root):
        self.groups[root][cage] = group
        self.filled_cages.add(cage)

    def __del__(self):
        logger = getLogger()
        logger.info("Completed.")


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
