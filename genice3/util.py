"""
Utility functions for genice3
"""

import numpy as np
from typing import List, Iterable
from logging import getLogger
from collections import defaultdict
import itertools as it

import networkx as nx
import pairlist as pl

from genice3.molecule import Molecule


def replicate_positions(positions1, replica_vectors, grand_cellmat):
    """レプリカ単位胞の数だけ、水分子位置を複製する。"""
    # 空の配列の場合は空の配列を返す
    if len(positions1) == 0:
        return np.array([]).reshape(0, 3)
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


def serialize(molecules: List[Molecule]):
    atoms = []
    for molecule in molecules:
        for name, position in zip(molecule.labels, molecule.sites):
            atoms.append([molecule.name, name, position])
    return atoms


# ============================================================================
# CIF-like data processing functions (moved from genice2.CIF)
# ============================================================================


def fullatoms(atomd, sops):
    """Generate all atoms from atom dictionary and symmetry operations."""
    global x, y, z
    logger = getLogger()
    # the variables to be evaluated by eval() must be global (?)
    atoms = []
    for sop in sops:
        for name, pos in atomd.items():
            x, y, z = pos
            p = sop(x, y, z)
            tooclose = False
            for n, f in atoms:
                d = f - p
                d -= np.floor(d + 0.5)
                L2 = d @ d
                if L2 < 0.0001:
                    logger.debug("Too close: {0} {1}".format(f, p))
                    tooclose = True
                    break
            if not tooclose:
                atoms.append((name, p))
    return atoms


def atomdic(atoms):
    """Parse atom positions from a string format."""
    atomd = dict()
    for atom in atoms.split("\n"):
        cols = atom.split()
        if len(cols) == 0:
            continue
        # remove error digits from 1-3 cols.
        xyz = []
        for col in cols[1:4]:
            p = col.find("(")
            if p >= 0:
                col = col[:p]
            xyz.append(float(col))
        name = cols[0]
        atomd[name] = xyz
    return atomd


def symmetry_operators(symops: str, offsets: Iterable = [("+0", "+0", "+0")]):
    """Generator of the symmetry operations

    Args:
        symops (str): The table of symmetry operators, e.g. http://img.chem.ucl.ac.uk/sgp/large/036az3.htm
        offset (Iterable, optional): Offset to the origin. Defaults to ["+0", "+0", "+0"].
    """

    def _wrap(x):
        return x
        return x - np.floor(x + 0.5)

    def eval_(v, x, y, z):
        return eval(v)

    symops = symops.translate(str.maketrans("XYZ", "xyz"))
    for symop in symops.split("\n"):
        cols = symop.split(",")
        if len(cols) <= 1:
            cols = symop.split()
            if len(cols) <= 1:
                continue
        for offset in offsets:
            ops = []
            for col, offs in zip(cols, offset):
                ops.append(col + offs)
            yield lambda x, y, z: _wrap(np.array([eval_(op, x, y, z) for op in ops]))


# def atoms_to_waters(oxygens, hydrogens, cell, partial_order=False):
#     logger = getLogger("atoms_to_waters")
#     oh = defaultdict(list)
#     parent = dict()
#     # find covalent OH bonds
#     for i, j in pl.pairs_iter(
#         oxygens, maxdist=0.15, cell=cell, pos2=hydrogens, distance=False  # nm
#     ):
#         oh[i].append(j)
#         parent[j] = i
#     logger.debug(parent)
#     pairs = []
#     # find HBs
#     for i, j in pl.pairs_iter(
#         oxygens, maxdist=0.20, cell=cell, pos2=hydrogens, distance=False  # nm
#     ):
#         if j not in oh[i]:
#             # H is on a different water molecule
#             p = parent[j]
#             pairs.append([p, i])

#     logger.debug(pairs)

#     if partial_order:
#         oo_pairs = [
#             (i, j)
#             for i, j in pl.pairs_iter(oxygens, maxdist=0.30, cell=cell, distance=False)
#         ]
#         # positions of the oxygen atoms, fixed edges, and unassigned edges.
#         return oxygens, pairs, oo_pairs

#     waters = []
#     for i in range(len(oh)):
#         logger.debug((i, oh[i]))
#         j, k = oh[i]
#         dhj = hydrogens[j] - oxygens[i]
#         dhk = hydrogens[k] - oxygens[i]
#         dhj -= np.floor(dhj + 0.5)
#         dhk -= np.floor(dhk + 0.5)
#         com = oxygens[i] + (dhj + dhk) / 18
#         waters.append(com)

#     return waters, pairs


def atoms_to_waters(oxygens, hydrogens, cell, partial_order=False):
    """
    原子座標 (O, H) から水分子の重心座標と水素結合ペアを求める。

    Parameters
    ----------
    oxygens : np.ndarray
        O原子の分数座標 (shape: (n_O, 3))
    hydrogens : np.ndarray
        H原子の分数座標 (shape: (n_H, 3))
    cell : np.ndarray
        セル行列
    partial_order : bool, default False
        True のときは、重心計算をせず、O–O ペアも返す（秩序化用）。

    Returns
    -------
    if partial_order is False:
        waters : list[np.ndarray]
            水分子の重心座標のリスト
        pairs : list[list[int, int]]
            水素結合ペア (donor_O_index, acceptor_O_index) のリスト

    if partial_order is True:
        oxygens : np.ndarray
        pairs   : list[list[int, int]]  (O–O 間のHB由来ペア)
        oo_pairs: list[tuple[int, int]] (距離しきい値以内の O–O ペア)
    """
    logger = getLogger("atoms_to_waters")

    # --- Step 1: 共有結合 OH（親Oを決める） ------------------------------------
    oh = defaultdict(list)  # Oインデックス -> 共有結合しているHインデックスのリスト
    parent = {}  # Hインデックス -> 親Oインデックス

    for o_idx, h_idx in pl.pairs_iter(
        oxygens, maxdist=0.15, cell=cell, pos2=hydrogens, distance=False  # nm
    ):
        oh[o_idx].append(h_idx)
        parent[h_idx] = o_idx

    logger.debug("covalent OH parent map: %s", parent)

    # --- Step 2: 水素結合ペア (O_donor, O_acceptor) を集める --------------------
    pairs = []
    for o_idx, h_idx in pl.pairs_iter(
        oxygens, maxdist=0.20, cell=cell, pos2=hydrogens, distance=False  # nm
    ):
        # 共有結合しているHは除外（別の水分子に属するHだけを使う）
        if h_idx not in oh[o_idx]:
            donor = parent[h_idx]  # このHの親O
            acceptor = o_idx
            pairs.append([donor, acceptor])

    logger.debug("HB pairs (donor, acceptor): %s", pairs)

    # --- Step 3: 秩序化用の部分情報だけ欲しい場合 ------------------------------
    if partial_order:
        oo_pairs = [
            (i, j)
            for i, j in pl.pairs_iter(oxygens, maxdist=0.30, cell=cell, distance=False)
        ]
        # oxygens: O原子座標
        # pairs  : HBペア (donor, acceptor)
        # oo_pairs: 近接 O–O ペア
        return oxygens, pairs, oo_pairs

    # --- Step 4: 各Oごとに水分子重心を計算 ------------------------------------

    def _water_com(o_index, h_indices):
        """1つのOと2つのHから、水分子重心の分数座標を計算する。"""
        h1_index, h2_index = h_indices

        o = oxygens[o_index]
        h1 = hydrogens[h1_index]
        h2 = hydrogens[h2_index]

        # 近接イメージを選ぶ（周期境界考慮）
        dh1 = h1 - o
        dh2 = h2 - o
        dh1 -= np.floor(dh1 + 0.5)
        dh2 -= np.floor(dh2 + 0.5)

        # O質量:16, H質量:1 として COM = (16*O + H1 + H2)/18
        return o + (dh1 + dh2) / 18.0

    waters = []
    for o_index in range(len(oh)):
        logger.debug("O %d has Hs %s", o_index, oh[o_index])
        # ここでは H がちょうど2個ある前提（既存コードと同じ前提）
        waters.append(_water_com(o_index, oh[o_index]))

    return waters, pairs


def waters_and_pairs(
    cell,
    atomd,
    sops,
    rep=(1, 1, 1),
    O_labels=("O",),
    H_labels="DH",
    partial_order=False,
):
    """Generate water positions and pairs from CIF-like data.

    Args:
        cell: Cell matrix
        atomd: Atom dictionary from atomdic()
        sops: Symmetry operators from symmetry_operators()
        rep: Replication factors (default: (1, 1, 1))
        O_labels: Labels for oxygen atoms (default: ("O",))
        H_labels: Labels for hydrogen atoms (default: "DH")
        partial_order: Whether to return partial order information (default: False)

    Returns:
        waters: Water positions (numpy array)
        pairs: Hydrogen bond pairs (list or None)
        oo_pairs: O-O pairs (only if partial_order=True)
    """
    # 部分秩序の場合は、位置が確定している水素だけをatomdに入れ、未確定の水素は省く。
    logger = getLogger()
    oxygens = []
    hydrogens = []
    for name, pos in fullatoms(atomd, sops):
        if name[0] in O_labels:
            oxygens.append(pos)
        elif name[0] in H_labels:
            hydrogens.append(pos)

    if not partial_order:
        assert (
            len(oxygens) * 2 == len(hydrogens) or len(hydrogens) == 0
        ), f"H {len(oxygens) * 2}: O {len(hydrogens)}"
    cell *= np.array(rep)

    oo = [
        [o[0] + x, o[1] + y, o[2] + z]
        for o in oxygens
        for x in range(rep[0])
        for y in range(rep[1])
        for z in range(rep[2])
    ]
    oxygens = np.array(oo)
    oxygens /= np.array(rep)

    # if no hydrogens are included in atomd,
    if len(hydrogens) == 0:
        return oxygens, None

    hh = [
        [h[0] + x, h[1] + y, h[2] + z]
        for h in hydrogens
        for x in range(rep[0])
        for y in range(rep[1])
        for z in range(rep[2])
    ]
    hydrogens = np.array(hh)
    hydrogens /= np.array(rep)

    return atoms_to_waters(oxygens, hydrogens, cell, partial_order=partial_order)


def shortest_distance(coord, cell):
    dmin = 1e99
    for i, j in it.combinations(coord, 2):
        d = i - j
        d -= np.floor(d + 0.5)
        d = d @ cell
        d = np.linalg.norm(d)
        if d < dmin:
            dmin = d
    return dmin


def density_in_g_cm3(nmol_water, cell):
    return 18 * nmol_water / (np.linalg.det(cell) * 1e-21 * 6.022e23)


def validate_ice_rules(dg: nx.DiGraph):
    logger = getLogger()
    valid = True
    for node in dg:
        in_degree = dg.in_degree(node)
        out_degree = dg.out_degree(node)
        if in_degree != 2 or out_degree != 2:
            valid = False
            logger.info(
                f"Node {node} has {dg.in_degree(node)} incoming edges"
                f" and {dg.out_degree(node)} outgoing edges"
            )
    return valid


# Candidates for ice XI
# {'P12_11', 'P2_12_12_1', 'Pca2_1', 'P1c1', 'Cmc2_1', 'C1c1', 'P12', 'Pbn2_1', 'Pna2_1'}
def operations(spaceg: str, origin: int = 0):
    symops = None
    if spaceg == "P12_11" or spaceg == "4":
        symops = symmetry_operators(
            """
            x,            y,            z
            -x,          1/2+y,         -z
        """
        )
    elif spaceg == "P2_12_12_1" or spaceg == "19":
        # http://img.chem.ucl.ac.uk/sgp/large/019a.htm
        if origin == 0:
            symops = symmetry_operators(
                """
                x,            y,            z
                1/2+x,        1/2-y,         -z
                -x,          1/2+y,        1/2-z
                1/2-x,         -y,          1/2+z
            """
            )
    elif spaceg == "Pca2_1" or spaceg == 29:
        symops = symmetry_operators(
            """
            x,            y,            z
            1/2-x,          y,          1/2+z
            1/2+x,         -y,            z
            -x,           -y,          1/2+z
        """
        )
    elif spaceg == "P1c1" or spaceg == "7":
        symops = symmetry_operators(
            """
            x,            y,            z
            x,           -y,          1/2+z
        """
        )
    elif spaceg == "Cmc2_1" or spaceg == "36":
        symops = symmetry_operators(
            """
            x,            y,            z
            -x,            y,            z
            x,           -y,          1/2+z
            -x,           -y,          1/2+z
        """,
            offsets=[("+0", "+0", "+0"), ("+1/2", "+1/2", "+0")],
        )
    elif spaceg == "C1c1" or spaceg == "9":
        symops = symmetry_operators(
            """
            x,            y,            z
            x,           -y,          1/2+z
        """,
            offsets=[("+0", "+0", "+0"), ("+1/2", "+1/2", "+0")],
        )
    elif spaceg == "P1" or spaceg == "1":
        symops = symmetry_operators(
            """
            x,            y,            z
        """
        )
    elif spaceg == "Pbn2_1":
        # there are two #33 (?)
        symops = symmetry_operators(
            """
            x,            y,            z
            1/2-x,        1/2+y,          z
            1/2+x,        1/2-y,        1/2+z
            -x,           -y,          1/2+z
        """
        )
    elif spaceg == "Pna2_1":
        # there are two #33 (?)
        symops = symmetry_operators(
            """
            x,            y,            z
            1/2-x,        1/2+y,        1/2+z
            1/2+x,        1/2-y,          z
            -x,           -y,          1/2+z
        """
        )
    return symops
