"""
Generate a hydrogen-disordered honeycomb bilayer ice.

Usage:
  genice2 bilayer                     Default size (6,6)
  genice2 bilayer[size=6,10]          Larger size
  genice2 bilayer[size=6,10:sw=0.2]   With Stone-Wales defects

Options:
  size=x,y   x must be a multiple of 3.
  sw=0.0     Z specifies the ratio of Stone-Wales defects to the number of
               sites.
"""

from genice2.cell import cellvectors
import genice2.lattices

import networkx as nx
import numpy as np

# import random # not a good manner to initialize its own random generator.
import itertools as it
from cycless.cycles import cycles_iter
from logging import getLogger


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {
    "ref": {"bilayer": "Koga 1997"},
    "usage": usage(),
    "brief": "A Bilayer Honeycomb Ice Phase in Hydrophobic Nanopores.",
    "test": ({"args": {"size": "6,10"}}, {"args": {"size": "15,18", "sw": "0.2"}}),
}


# Unnecessary when pairs are given.
# bondlen = 2.0 * 3.0 / 8.0**0.5 * 1.1

# import matplotlib.pyplot as plt


def remove_node(cycle, v):
    cycle.remove(v)


def insert_node(cycle, v, pair):
    a, b = pair
    apos = cycle.index(a)
    if cycle[apos - 1] == b:
        cycle.insert(apos, v)
    else:
        cycle.insert(apos + 1, v)


def stone_wales(g, cycles):
    while True:
        i = np.random.randint(len(g))
        ineis = set(g[i])
        j = np.random.choice(list(ineis))
        jneis = set(g[j])

        owners = g.edges[i, j]["owners"]

        ineis -= {j}
        jneis -= {i}

        j1, j2 = list(jneis)
        i1, i2 = list(ineis)
        for o1 in owners:
            cycle = cycles[o1]
            assert i in cycle
            assert j in cycle
            if i1 not in cycle:
                continue
            if j2 in cycle:
                j1, j2 = j2, j1
            break

        # i1 and j1 is in the same side
        # cycle o1 includes i1 and j1
        if owners[0] == o1:
            o2 = owners[1]
        else:
            o2 = owners[0]

        owners = g.edges[i, i1]["owners"]
        if owners[0] == o1:
            o3 = owners[1]
        else:
            o3 = owners[0]

        owners = g.edges[j, j1]["owners"]
        if owners[0] == o1:
            o4 = owners[1]
        else:
            o4 = owners[0]
        if (
            len(cycles[o1]) < 6
            or len(cycles[o2]) < 6
            or len(cycles[o3]) > 6
            or len(cycles[o4]) > 6
        ):
            continue
        break
    g.add_edge(i1, j)
    g.edges[i1, j]["owners"] = [o1, o3]
    g.add_edge(i, j2)
    g.edges[i, j2]["owners"] = [o2, o4]
    g.edges[i, j]["owners"] = [o3, o4]

    g.remove_edge(i, i1)
    g.remove_edge(j, j2)
    # print(cycles[o1])
    # print(cycles[o2])
    # print(cycles[o3])
    # print(cycles[o4])

    remove_node(cycles[o1], i)
    remove_node(cycles[o2], j)
    insert_node(cycles[o3], j, (i, i1))
    insert_node(cycles[o4], i, (j, j2))
    # print(cycles[o1])
    # print(cycles[o2])
    # print(cycles[o3])
    # print(cycles[o4])


# relaxation


def relax(g, rpos, cell, cycles=None, k=1.0):
    for loop in range(100):
        forces = np.zeros_like(rpos)
        for i, j in g.edges():
            d = rpos[i] - rpos[j]
            d -= np.floor(d + 0.5)
            d *= cell
            L = np.linalg.norm(d)
            f = k * (L - 1.0)
            f *= d / L
            forces[i] -= f
            forces[j] += f
        if cycles is not None:
            for cycle in cycles:
                for i, j in zip(cycle, cycle[-2:] + cycle[:-2]):
                    d = rpos[i] - rpos[j]
                    d -= np.floor(d + 0.5)
                    d *= cell
                    L = np.linalg.norm(d)
                    f = k * (L - 3**0.5)
                    f *= d / L
                    forces[i] -= f
                    forces[j] += f
        rpos += forces / cell


# def draw(g, rpos, cell):
#     nodes = rpos * cell
#     plt.axis('equal')
#     plt.scatter(nodes[:,0], nodes[:,1])
#     for i,j in g.edges():
#         d = rpos[i] - rpos[j]
#         s = np.floor(d+0.5)
#         if np.allclose(s, 0):
#             seg = np.vstack([nodes[j], nodes[j]+d*cell])
#             plt.plot(seg[:,0], seg[:,1], "-")
#         else:
#             d -= s
#             seg = np.vstack([nodes[j], nodes[j]+d*cell])
#             plt.plot(seg[:,0], seg[:,1], "-")
#             seg = np.vstack([nodes[i]-d*cell, nodes[i]])
#             plt.plot(seg[:,0], seg[:,1], "-")


class Lattice(genice2.lattices.Lattice):
    def __init__(self, **kwargs):
        logger = getLogger()

        NX, NY = 6, 6
        sw = 0.0

        for k, v in kwargs.items():
            if k == "size":
                NX, NY = [int(x) for x in v.split(",")]
                assert NX % 3 == 0, "X must be a multiple of 3."
            elif k == "sw":
                sw = float(v)
            elif v is True:
                usage()
                # unlabeled option
                sys.exit(1)

        logger.info(f"Bilayer ice of size ({NX}x{NY})")

        y1 = 3**0.5 / 2
        cell = np.array([NX, y1 * NY])

        nodes1 = [
            (i, j * y1 * 2) for i in range(NX) for j in range(NY // 2) if i % 3 != 2
        ]
        nodes2 = [
            (i + 0.5, (j + 0.5) * y1 * 2)
            for i in range(NX)
            for j in range(NY // 2)
            if i % 3 != 0
        ]
        nodes = np.array(nodes1 + nodes2)

        rpos = nodes / cell

        g = nx.Graph()
        for i, j in it.combinations(range(nodes.shape[0]), 2):
            d = rpos[i] - rpos[j]
            d -= np.floor(d + 0.5)
            d *= cell
            if d @ d < 1.01**2:
                g.add_edge(i, j)

        cycles = [list(cycle) for cycle in cycles_iter(g, 6)]

        for edge in g.edges():
            g.edges[edge]["owners"] = []
        for i, cycle in enumerate(cycles):
            for j, k in zip(cycle, cycle[-1:] + cycle[:-1]):
                g.edges[j, k]["owners"].append(i)

        # fig = plt.figure(figsize=(10,10))
        # draw(g, rpos, cell)
        # plt.show()

        # introduce defects
        Nnode = len(rpos)
        Nshuffle = int(Nnode * sw)
        # Nshuffle = 10
        for i in range(Nshuffle):
            stone_wales(g, cycles)

        relax(g, rpos, cell, cycles, k=0.1)

        self.cell = cellvectors(a=cell[0], b=cell[1], c=10.0) * 0.276

        self.density = Nnode * 2 * 18 / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
        self.waters = np.zeros([Nnode * 2, 3])
        self.waters[:Nnode, 0:2] = rpos
        self.waters[Nnode:, 0:2] = rpos
        self.waters[Nnode:, 2] = 0.276 / self.cell[2, 2]
        # self.bondlen = 1.1*0.276

        self.pairs = []
        for i, j in g.edges():
            self.pairs.append([i, j])
            self.pairs.append([i + Nnode, j + Nnode])
        for i in range(Nnode):
            self.pairs.append([i, i + Nnode])
