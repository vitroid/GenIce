"""
Utility functions for genice3
"""

import numpy as np
from typing import Dict, Tuple, List
from dataclasses import dataclass
import string
from logging import getLogger

import networkx as nx
from cycless.cycles import centerOfMass, cycles_iter
from cycless.polyhed import cage_to_graph, polyhedra_iter
from graphstat import GraphStat
from math import pi, acos, degrees

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
