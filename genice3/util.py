"""
Utility functions for genice3
"""

import numpy as np
from typing import Dict, Tuple, List, Iterable
from dataclasses import dataclass
import string
from logging import getLogger
from collections import defaultdict
import itertools as it

import networkx as nx
import pairlist as pl
from cycless.cycles import centerOfMass, cycles_iter
from cycless.polyhed import cage_to_graph, polyhedra_iter
from graphstat import GraphStat
from math import pi, acos, degrees, sin, cos

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


def cellshape(cellmat):
    """
    From cell matrix to a, b, c, alpha, beta, and gamma.
    """
    a = np.linalg.norm(cellmat[0])
    b = np.linalg.norm(cellmat[1])
    c = np.linalg.norm(cellmat[2])
    alpha = degrees(acos((cellmat[1] @ cellmat[2]) / (b * c)))
    beta = degrees(acos((cellmat[2] @ cellmat[0]) / (c * a)))
    gamma = degrees(acos((cellmat[0] @ cellmat[1]) / (a * b)))
    return a, b, c, alpha, beta, gamma


def cellvectors(a, b, c, A=90, B=90, C=90):
    """
    Generate cell vectors from a,b,c and alpha, beta, gamma.
    """
    # probably same as six2nine in rigid.py
    logger = getLogger()
    A *= pi / 180
    B *= pi / 180
    C *= pi / 180
    sA, cA = sin(A), cos(A)
    sB, cB = sin(B), cos(B)
    sC, cC = sin(C), cos(C)
    ea = np.array([1.0, 0.0, 0.0])
    eb = np.array([cC, sC, 0])
    # ec.ea = ecx = cos(B)
    # ec.eb = ecx*ebx + ecy*eby = cos(A)
    ecx = cB
    ecy = (cA - ecx * eb[0]) / eb[1]
    ecz = (1 - ecx**2 - ecy**2) ** 0.5
    ec = np.array([ecx, ecy, ecz])
    return np.vstack([ea * a, eb * b, ec * c])


@dataclass
class CageSpec:
    label: str  # A12, etc.
    faces: str  # 5^12 6^2, etc.
    graph: nx.Graph  # labels of water constituting the cage

    def to_json_capable_data(self):
        return {"label": self.label, "faces": self.faces, "nodes": list(self.graph)}

    def __repr__(self) -> str:
        return (
            f"CageSpec(label={self.label!r}, "
            f"faces={self.faces!r}, "
            f"n_nodes={self.graph.number_of_nodes()})"
        )

    def __str__(self) -> str:
        return f"{self.label} ({self.faces}) {self.graph.number_of_nodes()} nodes"


@dataclass
class CageSpecs:
    specs: list[CageSpec]
    positions: np.ndarray  # in fractional coordinates

    def to_json_capable_data(self):
        data = []
        for position, specs in zip(self.positions, self.specs):
            data.append(
                {
                    "frac_pos": position.tolist(),
                    "specs": specs.to_json_capable_data(),
                }
            )
        return dict(enumerate(data))

    def __repr__(self) -> str:
        return (
            f"CageSpecs(n_cages={len(self.specs)}, "
            f"positions_shape={self.positions.shape})"
        )


def _assign_unused_label(basename, labels):
    enum = 0
    label = f"A{basename}"
    while label in labels:
        char = string.ascii_lowercase[enum]
        label = f"A{basename}{char}"
        enum += 1
    return label


def _make_cage_expression(ring_ids, ringlist):
    ringcount = [0 for i in range(9)]
    for ring in ring_ids:
        ringcount[len(ringlist[ring])] += 1
    index = []
    for i in range(9):
        if ringcount[i] > 0:
            index.append(f"{i}^{ringcount[i]}")
    index = " ".join(index)
    return index


def assess_cages(graph, node_pos):
    """Assess cages from  the graph topology.

    Args:
        graph (graph-like): HB network
        nodepos (np.Array): Positions of the nodes
    """
    logger = getLogger()

    # Prepare the list of rings
    # taking the positions in PBC into account.
    ringlist = [
        [int(x) for x in ring] for ring in cycles_iter(nx.Graph(graph), 8, pos=node_pos)
    ]

    # Positions of the centers of the rings.
    ringpos = np.array([centerOfMass(ringnodes, node_pos) for ringnodes in ringlist])

    MaxCageSize = 22
    positions = []
    cagetypes = []
    # data storage of the found cages
    db = GraphStat()
    labels = set()
    g_id2label = dict()

    # Detect cages and classify
    cages = [cage for cage in polyhedra_iter(ringlist, MaxCageSize)]
    positions = [centerOfMass(list(cage), ringpos) for cage in cages]

    cagespecs = []
    cage_positions = []
    for cage, position in zip(cages, positions):
        g = cage_to_graph(cage, ringlist)
        cagesize = len(cage)
        g_id = db.query_id(g)
        # if it is a new cage type
        if g_id < 0:
            # new type!
            # register the last query
            g_id = db.register()

            # prepare a new label
            label = _assign_unused_label(cagesize, labels)
            g_id2label[g_id] = label
            labels.add(label)

            # cage expression
        else:
            label = g_id2label[g_id]
        faces = _make_cage_expression(cage, ringlist)
        cagespecs.append(CageSpec(label=label, faces=faces, graph=g))
        # print(f"{label=}, {faces=}, {ringlist=}")
        # print([len(ringlist[ring]) for ring in cage])

        cage_positions.append(position)
    if len(cage_positions) == 0:
        logger.info("    No cages detected.")
    return CageSpecs(specs=cagespecs, positions=np.array(cage_positions))


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

    oh = defaultdict(list)
    parent = dict()
    # find covalent OH bonds
    for i, j in pl.pairs_iter(
        oxygens, maxdist=0.15, cell=cell, pos2=hydrogens, distance=False  # nm
    ):
        oh[i].append(j)
        parent[j] = i
    logger.debug(parent)
    pairs = []
    # find HBs
    for i, j in pl.pairs_iter(
        oxygens, maxdist=0.20, cell=cell, pos2=hydrogens, distance=False  # nm
    ):
        if j not in oh[i]:
            # H is on a different water molecule
            p = parent[j]
            pairs.append([p, i])

    logger.debug(pairs)

    if partial_order:
        oo_pairs = [
            (i, j)
            for i, j in pl.pairs_iter(oxygens, maxdist=0.30, cell=cell, distance=False)
        ]
        # positions of the oxygen atoms, fixed edges, and unassigned edges.
        return oxygens, pairs, oo_pairs

    waters = []
    for i in range(len(oh)):
        logger.debug((i, oh[i]))
        j, k = oh[i]
        dhj = hydrogens[j] - oxygens[i]
        dhk = hydrogens[k] - oxygens[i]
        dhj -= np.floor(dhj + 0.5)
        dhk -= np.floor(dhk + 0.5)
        com = oxygens[i] + (dhj + dhk) / 18
        waters.append(com)

    return waters, pairs


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
