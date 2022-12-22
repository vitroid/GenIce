#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Graph variants with geometrical information
"""


from __future__ import print_function

import heapq
import math
import random
import sys
from logging import DEBUG, basicConfig, getLogger

import networkx as nx
import numpy as np
import yaplotlib as yp

# Dijkstra




def flatten(L):       # Flatten linked list of form [0,[1,[2,[]]]]
    while len(L) > 0:
        yield L[0]
        L = L[1]

# modified from http://code.activestate.com/recipes/119466/
# for networkx
# Extended for multiple goals


def shortest_path(G, start, ends):
    """
    Find a shortest path from the start to one of the ends.

    Returns:
    list of nodes on the shortest path from the start to one of the ends.
    """
    q = [(0, start, ())]  # Heap of (cost, path_head, path_rest).
    visited = set()       # Visited vertices.
    while True:
        try:
            ppap = heapq.heappop(q)
        except IndexError:
            return None
        (cost, v1, path) = ppap
        if v1 not in visited:
            visited.add(v1)
            if v1 in ends:
                return list(flatten(path))[::-1] + [v1]
            path = (v1, path)
            for v2 in G[v1]:
                if v2 not in visited and not G[v1][v2]['fixed']:
                    heapq.heappush(q, (cost + 1, v2, path))


def shortest_paths(G, start, ends, allowfixed=False):
    """
    Find all shortests paths from the start to one of the ends.

    Returns:
    list of shortest paths from the start to one of the ends.
    """
    # logger = getLogger()
    q = [(0, [start, ])]  # Heap of (cost, path)
    visited = set()       # Visited vertices.
    cheapest = 999999
    paths = []
    while len(q) > 0:
        # logger.debug(q)
        (cost, path) = heapq.heappop(q)
        if cost > cheapest:
            break
        v0 = path[-1]
        if v0 in ends:
            cheapest = cost  # first arrived
            paths.append(path)
        else:
            if v0 in visited:
                continue
            visited.add(v0)
            for v1 in G[v0]:
                if v1 not in visited:
                    if allowfixed or not G[v0][v1]['fixed']:
                        heapq.heappush(q, (cost + 1, path + [v1]))
    return paths


class YaplotDraw(nx.DiGraph):
    def __init__(self, coord, cell, data=None):
        super().__init__(data)
        self.coord = coord  # fractional coord.
        self.cell = cell

    def draw_edge(self, i, j):
        ci = self.coord[i]  # 0..1
        cj = self.coord[j]
        d = cj - ci
        d -= np.floor(d + 0.5)
        xi = ci @ self.cell
        xj = (ci + d) @ self.cell
        if self.has_edge(i, j):
            return yp.Color(4) + yp.ArrowType(2) + yp.Arrow(xi, xj)
        elif self.has_edge(j, i):
            return yp.Color(5) + yp.ArrowType(2) + yp.Arrow(xj, xi)
        else:
            return yp.Color(0) + yp.Line(xi, xj)

    def draw_cell(self):
        s = yp.Color(2)
        ex = np.array([1., 0., 0.])
        ey = np.array([0., 1., 0.])
        ez = np.array([0., 0., 1.])
        x = ex @ self.cell
        y = ey @ self.cell
        z = ez @ self.cell
        zero = np.zeros_like(x)
        for vx in (zero, x):
            for vy in (zero, y):
                s += yp.Line(vx + vy, vx + vy + z)
        for vx in (zero, x):
            for vz in (zero, z):
                s += yp.Line(vx + vz, vx + y + vz)
        for vz in (zero, z):
            for vy in (zero, y):
                s += yp.Line(vy + vz, x + vy + vz)
        return s

    def draw_path(self, path):
        s = yp.Color(3)
        for i in range(len(path) - 1):
            j, k = path[i], path[i + 1]
            s += self.draw_edge(j, k)
        return s

    def draw_all(self):
        s = yp.Color(3)
        s += self.draw_cell()
        for i, j in self.edges():
            s += self.draw_edge(i, j)
        return s


class IceGraph(nx.DiGraph):
    def __init__(self, data=None):
        super(IceGraph, self).__init__(data)
        # set of nodes that ignore the ice rule.
        # is added automatically in cationize/anionize
        self.immutables = set()

    def cationize(self, which):
        invert = set()
        fix = set()
        for i, j, data in self.edges(data=True):
            if j == which:
                invert.add((i, j))
                fix.add((j, i))
            elif i == which:
                fix.add((i, j))
        for i, j in invert:
            self.invert_edge(i, j)
        for i, j in fix:
            self[i][j]['fixed'] = True
        self.immutables.add(which)

    def anionize(self, which):
        invert = set()
        fix = set()
        for i, j, data in self.edges(data=True):
            if i == which:
                invert.add((i, j))
                fix.add((j, i))
            elif j == which:
                fix.add((i, j))
        for i, j in invert:
            self.invert_edge(i, j)
        for i, j in fix:
            self[i][j]['fixed'] = True
        self.immutables.add(which)

    def invert_edge(self, from_, to_, forced=False):
        """
        also invert the attribute vector
        """
        assert self.has_edge(from_, to_)
        assert forced or not self[from_][to_]['fixed']
        fix = self[from_][to_]['fixed']
        self.remove_edge(from_, to_)
        self.add_edge(to_, from_, fixed=fix)

    def invert_path(self, path, forced=False):
        for i in range(1, len(path)):
            f = path[i - 1]
            t = path[i]
            self.invert_edge(f, t, forced)  # also invert the attribute vector

    def invert_cycle(self, path, forced=False):
        p = path + [path[0], ]
        if not self.has_edge(p[0], p[1]):
            p = list(reversed(p))
        self.invert_path(p, forced=forced)

    def _goahead(self, node, marks, order):
        while node not in marks:
            marks[node] = len(order)
            order.append(node)
            nei = self.neighbors(node)
            next = random.randint(0, 1)
            assert len(nei) == 2, "Dangling bond: {0} {1}".format(node, nei)
            node = nei[next]
        # a cyclic path is found.
        # trim the uncyclic part
        order = order[marks[node]:]
        # add the node at the last again
        order.append(node)
        return order

    def homodromiccycle(self):
        """
        Randomly select a homodromic cycle
        """
        marks = dict()
        order = []
        node = random.randint(0, self.number_of_nodes() - 1)
        return self._goahead(node, marks, order)

    def isZ4(self):
        """
        Reply whether all the vertices have four neighbors or not.
        """
        undir = self.to_undirected()
        for node in range(undir.number_of_nodes()):
            if len(list(undir.neighbors(node))) != 4:
                return False
        return True

    def isZ22(self):
        """
        Reply whether all the vertices have two incomings and two outgoings.
        """
        for node in self:
            if len(list(self.successors(node))) != 2 or len(
                    list(self.predecessors(node))) != 2:
                return False
        return True

    def purgedefects(self, defects):
        """
        Buch's algorithm.
        """
        d = defects[0]
        # logger = getLogger()
        # logger.debug(self.immutables)
        if d in self.immutables:
            defects.pop(0)
            return
        if self.degree(d) != 4:  # TSL
            # assert self.degree(
            #     d) < 4, "Degree {0} should be <4.".format(self.degree(d))
            # logger.warn("  Defect {0} {1} >>{2} <<{3}".format(d,self.degree(d),self.in_degree(d),self.out_degree(d)))
            if self.degree(d) > 4 and self.out_degree(d) == 2:
                # accept if out-degree is 2.
                defects.pop(0)
                return
            if self.in_degree(d) <= 2 and self.out_degree(
                    d) <= 2:  # acceptable
                defects.pop(0)
                return
        if self.in_degree(d) == 2 and self.out_degree(d) == 2:
            defects.pop(0)
            return
        if self.in_degree(d) > 2:
            nodes = list(self.predecessors(d))
            i = random.randint(0, len(nodes) - 1)
            node = nodes[i]
            if not self[node][d]['fixed']:
                self.invert_edge(node, d)
                defects.append(node)
        if self.out_degree(d) > 2:
            nodes = list(self.successors(d))
            i = random.randint(0, len(nodes) - 1)
            node = nodes[i]
            if not self[d][node]['fixed']:
                self.invert_edge(d, node)
                defects.append(node)

    def bernal_fowler_defects(self):
        """
        Reply the list of defective vertices.

        It also counts the "ignore_ice_rules" sites.
        """
        # logger = getLogger()
        defects = []
        for i in range(self.number_of_nodes()):
            if self.in_degree(i) != 2 or self.out_degree(i) != 2:
                defects.append(i)
        return defects

    def excess_in_defects(self):
        """
        Reply the list of defective vertices.
        """
        for i in range(self.number_of_nodes()):
            if self.in_degree(i) > 2:
                yield i

    def excess_out_defects(self):
        """
        Reply the list of defective vertices.
        """
        for i in range(self.number_of_nodes()):
            if self.out_degree(i) > 2:
                yield i

    def purge_ice_defects(self):
        logger = getLogger()
        # TSL
        # if not self.isZ4():
        #    logger.error("Some water molecules do not have four HBs.")
        #    sys.exit(1)
        defects = self.bernal_fowler_defects()
        target = 1
        while len(defects) > target:
            target *= 2
        while len(defects) > 0:
            self.purgedefects(defects)
            if len(defects) <= target:
                logger.info("  Defects remaining: {0}".format(len(defects)))
                target //= 2
        # TSL
        # assert set(defects) == self.immutables, "Some water molecules do not obey the ice rule. {0} {1}".format(defects, self.immutables)

    def is_homodromic(self, path):
        for i in range(1, len(path)):
            if not self.has_edge(path[i - 1], path[i]):
                return False
        return True

    def is_cyclic_homodromic(self, path):
        p = path + [path[0], ]
        return self.is_homodromic(p) or self.is_homodromic(list(reversed(p)))


class SpaceIceGraph(IceGraph):
    """
    Digraph with geometrical info
    """
    XAXIS = 1
    YAXIS = 2
    ZAXIS = 3

    def __init__(self, data=None, coord=None, pbc=True, immutables=set()):
        # Both of them are slow for a huge system.
        super(SpaceIceGraph, self).__init__(data)
        if coord is not None:
            self.add_vectors(coord, pbc)  # fractional coord
        self.immutables = immutables

    def add_vectors(self, coord, pbc=True):
        """
        add vector attributes to each edge
        """
        self.coord = coord.copy()
        self.pbc = pbc
        for i, j, k in self.edges(data=True):
            vec = coord[j] - coord[i]
            if pbc:
                vec -= np.floor(vec + 0.5)
            k["vector"] = vec  # each edge has "vector" attribute

    def digraph(self):
        """
        Return the digraph only.
        """
        di = nx.DiGraph()
        for i, j, k in self.edges(data=True):
            di.add_edge(i, j, **k)
        return di

    def dipole_moment(self, order):
        """
        normally zero for a cycle.
        Non-zero when the cycle goes across the cell.

        For a cycle, the first element of the order must be the same as the last one.
        """
        delta = 0
        for i in range(len(order) - 1):
            delta += self.get_edge_data(order[i], order[i + 1])["vector"]
        return delta

    def invert_edge(self, from_, to_, forced=False):
        """
        also invert the attribute vector
        """
        assert self.has_edge(from_, to_)
        v = self.get_edge_data(from_, to_)["vector"]
        fixed = self.get_edge_data(from_, to_)["fixed"]
        assert forced or not fixed
        self.remove_edge(from_, to_)
        self.add_edge(to_, from_, vector=-v, fixed=fixed)

    def net_polarization(self):
        dipole = 0
        for i, j, k in self.edges(data=True):
            dipole += k["vector"]
        return dipole

    def vector_check(self):
        logger = getLogger()
        for i, j, k in self.edges(data=True):
            if k is None:
                logger.error("The edge ({0},{1}) has no vector.".format(i, j))


def find_apsis(vertices, coord, cell, distance, vertex, axis):
    # logger = getLogger()
    # for Z case
    apsis = coord[vertex] + axis * 0.5
    # find the atoms near the apsis
    min_r = 1e99
    min_a = -1
    for i in vertices:
        p = coord[i]
        d = p - apsis
        d -= np.floor(d + 0.5)
        a = d @  cell
        r = a @ a
        if r < min_r:
            min_a = i
            min_r = r
    return min_a


def estimate_edge_length(spaceicegraph, cell, vertex):
    logger = getLogger()
    # In case an anion is selected by bad fortune.
    if len(spaceicegraph.adj[vertex]) == 0:
        return 0
    nei = list(spaceicegraph.adj[vertex])[0]
    delta = spaceicegraph.coord[vertex] - spaceicegraph.coord[nei]
    delta -= np.floor(delta + 0.5)
    delta = delta @ cell  # distance in abs coord
    distance = np.linalg.norm(delta)
    logger.debug("Distance={0}".format(distance))
    return distance


def traversing_cycles_iter(spaceicegraph, cell, axis):
    """
    Find a farthest atom from the given atom, and
    make the shortest paths between them.
    """
    logger = getLogger()

    Nnode = spaceicegraph.number_of_nodes()
    # 近傍の定義のために結合距離が必要。
    distance = 0
    while distance == 0:
        vertex = random.sample(spaceicegraph.nodes(), 1)[0]
        distance = estimate_edge_length(spaceicegraph, cell, vertex)
    # すべての頂点を順番にあたる。
    for vertex in random.sample(spaceicegraph.nodes(), Nnode):
        apsis = find_apsis(spaceicegraph.nodes(
        ), spaceicegraph.coord, cell, distance * 1.3, vertex, axis)
        logger.debug("Apsis of {0}: {1}".format(vertex, apsis))
        path1 = shortest_path(spaceicegraph, vertex, [apsis, ])
        logger.debug("Path1: {0}".format(path1))
        paths = shortest_paths(spaceicegraph, vertex, [apsis, ])
        logger.debug("paths1: {0}".format(paths))
        if path1 is None:
            # No path found, probably because of the double networks
            continue
        logger.debug("Dipole of the path1: {0}".format(
            spaceicegraph.dipole_moment(path1)))
        path2 = shortest_path(spaceicegraph, apsis, [vertex, ])
        logger.debug("Path2: {0}".format(path2))
        if path2 is None:
            # No path found, probably because of the double networks
            continue
        logger.debug("Dipole of the path2: {0}".format(
            spaceicegraph.dipole_moment(path2)))
        # they should make a cycle.
        # It should go across the cell with a probability of 0.5.
        # It should g across the cell in an expected direction
        # with a probability of 0.25.
        # So we need some loops to get the expected one.
        cycle = path1 + path2[1:]
        d = spaceicegraph.dipole_moment(cycle) - axis
        logger.debug("Axis: {0} {1}".format(
            axis, spaceicegraph.dipole_moment(cycle)))
        rr = d @ d
        if rr < 0.1:
            logger.debug("Dipole of the harvest: {0}".format(
                spaceicegraph.dipole_moment(cycle)))
            yield cycle


def depolarize_(
        subgraph,
        coord,
        immutables=[],
        pbc=True,
        cell=np.identity(3),
        depol="strict"):
    """
    Receives a connected component of a large graph.
    Find a farthest atom (apsis) from the given atom, and
    make the shortest paths between them.

    [depol]
        strict:  die if pol is non-zero
        optimal: make pol smallest
        none:    do nothing

    """
    logger = getLogger()

    if depol == "none":
        logger.info("  Skip depolarization by request.")
        return subgraph

    spaceicegraph = SpaceIceGraph(
        subgraph,
        coord=coord,
        pbc=pbc,
        immutables=immutables)
    spaceicegraph.vector_check()

    # TSL
    # defect-defect chains
    defects = []
    for node in spaceicegraph.nodes():
        if spaceicegraph.degree(node) != 4:
            # assert spaceicegraph.degree(node) < 4
            defects.append(node)
    logger.info(f"  Water molecules that do not have 4 neighbors: {defects}")
    if len(defects) > 0:
        reject_count = 10
        while reject_count > 0:
            net_polar = spaceicegraph.net_polarization()
            if net_polar @ net_polar < 0.05**2:
                break  # without problem
            while True:
                orig = defects[random.randint(0, len(defects) - 1)]
                if spaceicegraph.out_degree(orig) == 2:
                    break
            while True:
                dest = defects[random.randint(0, len(defects) - 1)]
                if spaceicegraph.in_degree(dest) == 2:
                    break
            path = shortest_path(spaceicegraph, orig, [dest, ])
            if path is None:
                continue
            logger.debug("  Dipole moment = {0}".format(
                spaceicegraph.dipole_moment(path)))
            new_net_polar = net_polar - spaceicegraph.dipole_moment(path) * 2
            if np.linalg.norm(new_net_polar) < np.linalg.norm(net_polar):
                spaceicegraph.invert_path(path)
                net_polar = new_net_polar
                logger.info(
                    "  Net polarization: [{0:.2f} {1:.2f} {2:.2f}]".format(
                        *net_polar))
                reject_count = 10
            else:
                logger.debug("  Reject inversion")
                reject_count -= 1

    for axis_ in ([1, 0, 0], [-1, 0, 0], [0, 1, 0],
                  [0, -1, 0], [0, 0, 1], [0, 0, -1]):
        net_polar = spaceicegraph.net_polarization()
        logger.debug(f"Net polarization: {net_polar}")
        axis = np.array(axis_, dtype=float)
        L1 = (net_polar**2).sum()
        L2 = ((net_polar - axis * 2)**2).sum()
        if L2 < L1:
            logger.info(
                "  Net polarization: [{0:.2f} {1:.2f} {2:.2f}]".format(
                    *net_polar))
            for cycle in traversing_cycles_iter(spaceicegraph, cell, axis):
                if cycle is not None:
                    edges = [(cycle[i], cycle[i + 1])
                             for i in range(len(cycle) - 1)]
                    if len(edges) != len(set(edges)):
                        logger.debug("The cycle is entangled.")
                    else:
                        spaceicegraph.invert_path(cycle)
                        spaceicegraph.vector_check()
                        net_polar = spaceicegraph.net_polarization()
                        logger.info(
                            "  Net polarization: [{0:.2f} {1:.2f} {2:.2f}]".format(
                                *net_polar))
                        L1 = (net_polar**2).sum()
                        L2 = ((net_polar - axis * 2)**2).sum()
                        if L2 > L1:
                            break
            # break this for-loop

    if depol == "strict":
        if not np.allclose(net_polar, 0):
            logger.error(
                "  Gave up on total depolarization. Perhaps because they contain ions"
                " or the cells are very small."
                " Use --depol=optimal or --depol=none instead.")
            sys.exit(1)
    return spaceicegraph.digraph()


def depolarize(
        digraph,
        coord,
        immutables=[],
        pbc=True,
        cell=np.identity(3),
        depol="strict"):
    logger = getLogger()
    graph = nx.Graph(digraph)
    newdigraph = nx.DiGraph()
    for i, subgraph in enumerate(nx.connected_components(graph)):
        logger.info(f"  Component {i}:")
        subdigraph = digraph.subgraph(subgraph)  # really?
        newsubdigraph = depolarize_(
            subdigraph,
            coord=coord,
            immutables=immutables,
            pbc=pbc,
            cell=cell,
            depol=depol)
        newdigraph = nx.union(newdigraph, newsubdigraph)
    return newdigraph


def purge_ice_defects(icegraph):
    """
    This is faster than the method in icegraph, but
    it also polarizes the graph in the course of purging.
    """
    logger = getLogger()
    while len(icegraph.bernal_fowler_defects()) > 0:
        logger.info("# of defects: {0}".format(
            len(icegraph.bernal_fowler_defects())))
        ins = set(icegraph.excess_in_defects())
        for out in icegraph.excess_out_defects():
            while icegraph.out_degree(out) > 2:
                path = shortest_path(icegraph, out, [ins, ])
                if path is not None:
                    logger.debug("# of in defects: {0}".format(len(ins)))
                    icegraph.invert_path(path)
                    end = path[-1]
                    if icegraph.in_degree(end) == 2:
                        ins.remove(end)
                # logger.debug("IN:{0}".format(len(set(icegraph.excess_in_defects()))))
                # logger.debug("OUT:{0}".format(len(set(icegraph.excess_out_defects()))))


def test():
    # logger
    basicConfig(
        level=DEBUG,
        format='%(asctime)s- %(name)s - %(levelname)s - %(message)s')
    logger = getLogger(__name__)
    logger.setLevel(DEBUG)
    g = nx.Graph()
    # 6-cycle
    for i in range(5):
        g.add_edge(i, i + 1)
    g.add_edge(0, 5)
    print(shortest_paths(g, 0, [3, ], allowfixed=True))


if __name__ == "__main__":
    test()
