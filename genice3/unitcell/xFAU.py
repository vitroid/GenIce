#!/usr/bin/env python3
"""
Generate a (hydrogen-ordered) Aeroice nxFAU.

Usage:
  genice xFAU[3] > 3xfau.gro

Options:
  n    Length of the hexagonal prism. (rep=1:FAU, rep=0:SOD, rep>1: aeroice)
"""


from genice3.unitcell import ice1c as ic  # base topology
import genice3.unitcell
from cif2ice import cellvectors
from genice3.util import shortest_distance
import re
from logging import getLogger
import numpy as np
from collections import defaultdict
from math import acos, pi, sin, cos
import networkx as nx
import pairlist as pl

desc = {
    "ref": {"xFAU": "Matsui 2017"},
    "usage": __doc__,
    "brief": "Aeroice xFAU.",
    "test": ( # genice3-style options 
        {"rep": "0"},
        {"rep": "1"},
        {"rep": "2"},
        {"rep": "4"},
        {"rep": "8"},
    ),
}
# FAU Decoration of a 4-network
# 読みこんだAR3Rの座標を、FAU構造における多面体vertexの位置とみなし、
# それらを連結するネットワークを六角柱で修飾して大きなネットワークを作る。


# Aeroiceの超格子の接点多面体で、きちんとorderするように設計する。
# そのためには、超格子のもととなるdiamond latticeを白黒二部グラフとし、
# 黒から白へ向けていつも多角柱を作るようにする。
# 黒、白の接点多面体の六員環はすべてhomodromicとするが、どっちむきになるかは角柱の長さによって変わる。角柱が奇数段であれば、白黒とも同じ向きになるが、偶数段だと反転する。tune_anglesはてきとうに角柱を回転してつじつまをあわしているが、これを廃止する必要があるな。まずは全部決定論的に作る。
# ちょっとまじめに考えないといけない感じ。模型を睨む。


def tune_angles(sixvecs, pivot):
    """
    Find the best origin of angles to make sum cos(3 th) largest
    """
    sixangles = []
    for i in range(len(sixvecs)):
        vec = sixvecs[i]
        cosine = sixvecs[0] @ vec
        if cosine > 1.0:
            cosine = 1.0
        angle = acos(cosine)
        sine = np.cross(sixvecs[0], vec)
        if sine @ pivot < 0:
            angle = -angle
        sixangles.append(angle)
    offset = 0
    while True:
        sum = 0.0
        dsum = 0.0
        for a in sixangles:
            sum += cos((a + offset) * 6)
            dsum += -sin((a + offset) * 6)
        doffset = dsum / 20.0
        if abs(doffset) < 1e-6:
            return offset
        offset += doffset


class decorate:
    def __init__(self, atoms, cell, graph, Ncyl):
        """
        Ncyl is the number of cylinders to be inserted (>=0)
        """
        # make netghbor list
        self.diamond_graph = graph
        self.diamond_vertices = atoms
        self.diamond_cell = cell
        self.Ncyl = Ncyl
        self.vertices = []
        self.fixedEdges = []
        for edge in graph.edges():
            self.hexagonal_piller(edge)

    def hexagonal_piller(self, edge):
        logger = getLogger()
        i, j = edge
        dij = self.diamond_vertices[j] - self.diamond_vertices[i]
        dij -= np.floor(dij + 0.5)
        dij = dij @ self.diamond_cell
        scale = np.linalg.norm(dij)
        # dijは六角軸方向の単位ベクター
        dij /= scale
        # restsはiの隣接ノードのうちjを除いたもの
        rests = list(self.diamond_graph.neighbors(i))
        rests.remove(j)
        logger.debug(f"Rests: {rests}")
        # print(f"Rests: {list(rests)}")
        # Regularize the dihedral angles
        # to point them 6-fold directions.
        # by adding an offset
        vecs = []
        for k in rests:
            vec = self.diamond_vertices[k] - self.diamond_vertices[i]
            vec -= np.floor(vec + 0.5)
            vec = vec @ self.diamond_cell
            # orthogonalize
            shadow = dij @ vec
            vec -= shadow * dij
            vec /= np.linalg.norm(vec)
            vecs.append(vec)
        # print(f"{vecs=}")
        # 向きを同じにする。
        if np.linalg.det(np.vstack([dij, vecs[0], vecs[1]])) < 0:
            vecs[0], vecs[1] = vecs[1], vecs[0]
        offset = pi / 6  # 30 degree
        x = vecs[0]
        z = dij
        y = np.cross(z, x)
        sixvecs = np.zeros((6, 3))
        for j in range(6):
            a = j * pi * 2 / 6 + offset
            sixvecs[j] = x * cos(a) + y * sin(a)
        # determine r
        # assume edge length is 1
        # the radius of the outer sphere of the polyhed is sqrt(3/2)
        L = (3 / 2) ** 0.5 * 2 + self.Ncyl
        r = 1 / L  # edge len = radius of cyl
        rp = (3 / 2) ** 0.5 / L  # = radius of polyhed
        #
        icell = np.linalg.inv(self.diamond_cell)
        a = self.diamond_vertices[i] @ self.diamond_cell
        s = ""
        for j in range(0, self.Ncyl + 1):
            vec0 = dij * (rp + j * r) * scale + a
            for vec in sixvecs:
                rpos = vec0 + vec * r * scale
                pos = rpos @ icell
                # 水分子の位置。
                self.vertices.append(pos)

            first = len(self.vertices) - 6
            # 六角の結合を追加する。
            if j % 2 == 0:
                for k in range(5):
                    self.fixedEdges.append((first + k, first + k + 1))
                self.fixedEdges.append((first + 5, first))
            else:
                for k in range(5):
                    self.fixedEdges.append((first + k + 1, first + k))
                self.fixedEdges.append((first, first + 5))
            # 縦方向の結合を追加する交互に。
            if j > 0:
                for k in range(6):
                    if k % 2 == 0:
                        self.fixedEdges.append((first + k, first + k - 6))
                    else:
                        self.fixedEdges.append((first + k - 6, first + k))


class UnitCell(genice3.unitcell.UnitCell):
    """
    Generate a (hydrogen-ordered) Aeroice nxFAU.

    Options:
      rep=n    Length of the hexagonal prism. (rep=1:FAU, rep=0:SOD, rep>1: aeroice)
    """

    def __init__(self, **kwargs):
        logger = getLogger()
        ice1c = ic.UnitCell()
        cell1c = ice1c.cell
        waters1c = ice1c.lattice_sites
        graph1c = ice1c.graph
        #
        # 0..3を黒、4..7を白とする。もともと二部グラフになっているようだ。
        #
        Ncyl = None
        if len(kwargs) == 1:
            if "rep" in kwargs:
                Ncyl = int(kwargs["rep"])
            else:
                for k, v in kwargs.items():
                    if re.match("^[0-9]+$", k) is not None and v is True:
                        Ncyl = int(k)
                        break
        if Ncyl is None:
            raise ValueError("rep=n is required")
        logger.info("Superlattice {0}xFAU".format(Ncyl))
        dec = decorate(waters1c, cell1c, graph1c, Ncyl)

        coord = "relative"
        cell = cellvectors(a=dec.diamond_cell[0, 0], b=dec.diamond_cell[1, 1], c=dec.diamond_cell[2, 2])
        lattice_sites = np.array(dec.vertices)
        lattice_sites -= np.floor(lattice_sites)

        # 自力で密度を推定する。
        bondlen = shortest_distance(lattice_sites, cell)
        # 結合の長さはどんなに短くても0.276 nm。
        cell *= 0.276 / bondlen
        
        graph = nx.Graph([(i, j) for i, j in pl.pairs_iter(lattice_sites, maxdist=0.3, cell=cell, distance=False)])

        logger.info("after pairs_iter")
        super().__init__(
            cell=cell,
            lattice_sites=lattice_sites,
            coord=coord,
            fixed=nx.DiGraph(dec.fixedEdges),
            graph=graph,
        )
