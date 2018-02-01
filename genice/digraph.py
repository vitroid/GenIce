#!/usr/bin/env python
#coding: utf-8
#input: coordinate of the nodes, the digraph obeying the ice rule.
#output: the digraph with zero net dipole.

from __future__ import print_function
import sys
import math
import networkx
import random
import numpy as np
import logging
from genice import yaplotlib as yp
from genice import pairlist as pl

#Dijkstra


import heapq

def flatten(L):       # Flatten linked list of form [0,[1,[2,[]]]]
    while len(L) > 0:
        yield L[0]
        L = L[1]

#modified from http://code.activestate.com/recipes/119466/
#for networkx
#Extended for multiple goals
def shortest_path(G, start, ends):
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





class YaplotDraw(networkx.DiGraph):
    def __init__(self, coord, cell, data=None):
        super().__init__(data)
        self.coord = coord # fractional coord.
        self.cell  = cell

    def draw_edge(self,i,j):
        ci = self.coord[i]  #0..1
        cj = self.coord[j]
        d = cj - ci
        d -= np.floor( d + 0.5 )
        xi = np.dot(ci,self.cell)
        xj = np.dot(ci+d,self.cell)
        if self.has_edge(i,j):
            return yp.Color(4) + yp.ArrowType(2) + yp.Arrow(xi, xj)
        elif self.has_edge(j,i):
            return yp.Color(5) + yp.ArrowType(2) + yp.Arrow(xj, xi)
        else:
            return yp.Color(0) + yp.Line(xi, xj)
            

    def draw_cell(self):
        s = yp.Color(2)
        ex = np.array([1.,0.,0.])
        ey = np.array([0.,1.,0.])
        ez = np.array([0.,0.,1.])
        x = np.dot(ex,self.cell)
        y = np.dot(ey,self.cell)
        z = np.dot(ez,self.cell)
        zero = np.zeros_like(x)
        for vx in (zero, x):
            for vy in (zero, y):
                s += yp.Line(vx+vy,vx+vy+z)
        for vx in (zero, x):
            for vz in (zero, z):
                s += yp.Line(vx+vz,vx+y+vz)
        for vz in (zero, z):
            for vy in (zero, y):
                s += yp.Line(vy+vz,x+vy+vz)
        return s

    def draw_path(self,path):
        s = yp.Color(3)
        for i in range(len(path)-1):
            j,k = path[i],path[i+1]
            s += self.draw_edge(j,k)
        return s
    

    def draw_all(self):
        s = yp.Color(3)
        s += self.draw_cell()
        for i,j in self.edges():
            s += self.draw_edge(i,j)
        return s


class IceGraph(networkx.DiGraph):
    def __init__(self, data=None):
        super(IceGraph, self).__init__(data)
        #set of nodes that ignore the ice rule.
        #is added automatically in cationize/anionize
        self.ignores = set()  

#    def register_pairs(self,pairs):
#        self.clear()
#        for pair in pairs:
#            x,y = pair[:2]
#            self.add_edge(x,y)

    def cationize(self, which):
        invert = set()
        fix    = set()
        for i,j,data in self.edges(data=True):
            if j==which:
                invert.add((i,j))
                fix.add((j,i))
            elif i==which:
                fix.add((i,j))
        for i,j in invert:
            self.invert_edge(i,j)
        for i,j in fix:
            self[i][j]['fixed'] = True
        self.ignores.add(which)

    def anionize(self, which):
        invert = set()
        fix    = set()
        for i,j,data in self.edges(data=True):
            if i==which:
                invert.add((i,j))
                fix.add((j,i))
            elif j==which:
                fix.add((i,j))
        for i,j in invert:
            self.invert_edge(i,j)
        for i,j in fix:
            self[i][j]['fixed'] = True
        self.ignores.add(which)
        
    
    def invert_edge(self,from_,to_):
        """
        also invert the attribute vector
        """
        if not self.has_edge(from_, to_):
            logging.getLogger().error("No edge ({0},{1}).".format(from_,to_))
        assert not self[from_][to_]['fixed']
        self.remove_edge(from_,to_)
        self.add_edge(to_,from_,fixed=False)

        
    def invert_path(self, path):
        for i in range(len(path)-1):
            f = path[i]
            t = path[i+1]
            self.invert_edge(f,t) #also invert the attribute vector
            
    def _goahead(self,node,marks,order):
        while node not in marks:
            marks[node] = len(order)
            order.append(node)
            nei = self.neighbors(node)
            next = random.randint(0,1)
            assert len(nei) == 2, "Dangling bond: {0} {1}".format(node,nei)
            node = nei[next]
        #a cyclic path is found.
        #trim the uncyclic part
        order = order[marks[node]:]
        #add the node at the last again
        order.append(node)
        return order


    def homodromiccycle(self):
        """
        Randomly select a homodromic cycle
        """
        marks = dict()
        order = []
        node = random.randint(0,self.number_of_nodes()-1)
        return self._goahead(node,marks,order)


    def isZ4(self):
        """
        Reply whether all the vertices have four neighbors or not.
        """
        undir = self.to_undirected()
        for node in range(undir.number_of_nodes()):
            if len(list(undir.neighbors(node))) != 4:
                return False
        return True


    def purgedefects(self, defects):
        d = defects[0]
        logger = logging.getLogger()
        logger.debug(self.ignores)
        if d in self.ignores:
            defects.pop(0)
            return
        if self.degree(d) != 4:  # TSL
            assert self.degree(d) < 4
            # logger.warn("  Defect {0} {1} >>{2} <<{3}".format(d,self.degree(d),self.in_degree(d),self.out_degree(d)))
            if self.in_degree(d) <= 2 and self.out_degree(d) <= 2:  #acceptable
                defects.pop(0)
                return
        if self.in_degree(d) == 2 and self.out_degree(d) == 2:
            defects.pop(0)
            return
        if self.in_degree(d) > 2:
            nodes = list(self.predecessors(d))
            i = random.randint(0,len(nodes)-1)
            node = nodes[i]
            if not self[node][d]['fixed']:
                self.invert_edge(node,d)
                defects.append(node)
        if self.out_degree(d) > 2:
            nodes = list(self.successors(d))
            i = random.randint(0,len(nodes)-1)
            node = nodes[i]
            if not self[d][node]['fixed']:
                self.invert_edge(d,node)
                defects.append(node)

            
    def bernal_fowler_defects(self):
        """
        Reply the list of defective vertices.

        It also counts the "ignore_ice_rules" sites.
        """
        logger = logging.getLogger()
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
        logger = logging.getLogger()
        # TSL
        # if not self.isZ4():
        #    logger.error("Some water molecules do not have four HBs.")
        #    sys.exit(1)
        defects = self.bernal_fowler_defects()
        target = 1
        while len(defects) > target:
            target *= 2
        while len(defects)>0:
            self.purgedefects(defects)
            if len(defects) <= target:
                logger.info("  Defects remaining: {0}".format(len(defects)))
                target //= 2
        # TSL
        # assert set(defects) == self.ignores, "Some water molecules do not obey the ice rule. {0} {1}".format(defects, self.ignores)

                
    def is_homodromic(self, path):
        for i in range(len(path)-1):
            if not self.has_edge(path[i], path[i+1]):
                return False
        return True
            

            


class SpaceIceGraph(IceGraph):
    """
    Digraph with geometrical info
    """
    XAXIS=1
    YAXIS=2
    ZAXIS=3
    def __init__(self, data=None, coord=None, pbc=True, ignores=set()):
        super(SpaceIceGraph, self).__init__(data)
        if coord is not None:
            self.add_vectors(coord, pbc) # fractional coord
        self.ignores = ignores
            
    def add_vectors(self, coord, pbc=True):
        """
        add vector attributes to each edge
        """
        self.coord = coord #Shall it be copied?
        for i,j,k in self.edges(data=True):
            vec = coord[j] - coord[i]
            if pbc:
                vec -= np.floor(vec + 0.5)
            k["vector"] = vec  #each edge has "vector" attribute
        
    def dipole_moment(self, order):
        """
        normally zero for a cycle.
        Non-zero when the cycle goes across the cell.

        For a cycle, the first element of the order must be the same as the last one.
        """
        delta = np.zeros(3)
        for i in range(len(order)-1):
            delta += self.get_edge_data(order[i],order[i+1])["vector"]
        return delta
            
    def invert_edge(self,from_,to_):
        """
        also invert the attribute vector
        """
        if not self.has_edge(from_, to_):
            logging.getLogger().error("No edge ({0},{1}).".format(from_,to_))
        v = self.get_edge_data(from_,to_)["vector"]
        fixed = self.get_edge_data(from_,to_)["fixed"]
        assert not fixed
        self.remove_edge(from_,to_)
        self.add_edge(to_,from_,vector=-v, fixed=False)


        
    def net_polarization(self):
        dipole = np.zeros(3)
        for i,j,k in self.edges(data=True):
            dipole += k["vector"]
        return dipole

    def vector_check(self):
        logger = logging.getLogger()
        for i,j,k in self.edges(data=True):
            if k is None:
                logger.error("The edge ({0},{1}) has no vector.".format(i,j))
                


def find_apsis(coord, cell, distance, vertex, axis):
    logger = logging.getLogger()
    grid = pl.determine_grid(cell, distance)
    logger.debug("Grid: {0}".format(grid))
    #for Z case
    apsis = coord[vertex] + axis*0.5
    #find the atoms near the apsis
    pairs = pl.pairlist_fine_hetero([apsis,], coord, distance, cell, grid, distance=True)
    logger.debug("Neighbors of the apsis: {0}".format(pairs))
    min_d = 1e99
    min_a = -1
    for i,j,d in pairs:
        if d < min_d:
            min_a = j
            min_d = d
    return min_a
        

def estimate_edge_length(spaceicegraph, cell, vertex):
    logger = logging.getLogger()
    # In case an anion is selected by bad fortune.
    if len(spaceicegraph.adj[vertex]) == 0:
        return 0
    nei = list(spaceicegraph.adj[vertex])[0]
    delta = spaceicegraph.coord[vertex] - spaceicegraph.coord[nei]
    delta -= np.floor( delta + 0.5 )
    delta = np.dot(delta, cell)  #distance in abs coord
    distance = np.dot(delta,delta)**0.5
    logger.debug("Distance={0}".format(distance))
    return distance


def traversing_cycle(spaceicegraph, cell, axis, draw=None):
    """
    Find a farthest atom from the given atom, and
    make the shortest paths between them.
    """
    logger = logging.getLogger()

    distance = 0
    while distance == 0:
        vertex = random.randint(0,spaceicegraph.number_of_nodes()-1)
        distance = estimate_edge_length(spaceicegraph, cell, vertex)
    while True:
        vertex = random.randint(0,spaceicegraph.number_of_nodes()-1)
        apsis = find_apsis(spaceicegraph.coord, cell, distance*1.3, vertex, axis)
        logger.debug("Apsis of {0}: {1}".format(vertex, apsis))
        path1 = shortest_path(spaceicegraph, vertex, [apsis,])
        logger.debug("Path1: {0}".format(path1))
        if path1 is None:
            #No path found, probably because of the double networks
            continue
        logger.debug("Dipole of the path1: {0}".format(spaceicegraph.dipole_moment(path1)))
        path2 = shortest_path(spaceicegraph, apsis, [vertex,])
        logger.debug("Path2: {0}".format(path2))
        if path2 is None:
            #No path found, probably because of the double networks
            continue
        logger.debug("Dipole of the path2: {0}".format(spaceicegraph.dipole_moment(path2)))
        #they should make a cycle.
        #It should go across the cell with a probability of 0.5.
        #It should g across the cell in an expected direction
        #with a probability of 0.25.
        #So we need some loops to get the expected one.
        cycle = path1 + path2[1:]
        d = spaceicegraph.dipole_moment(cycle) - axis
        logger.debug("Axis: {0} {1}".format(axis,spaceicegraph.dipole_moment(cycle)))
        rr = np.dot(d,d)
        if rr < 0.1:
            break
    logger.debug("Dipole of the harvest: {0}".format(spaceicegraph.dipole_moment(cycle)))
    return cycle

def depolarize(spaceicegraph, cell, draw=None):
    """
    Find a farthest atom (apsis) from the given atom, and
    make the shortest paths between them.

    It works much better than depolarize()
    """
    logger = logging.getLogger()
    logger.debug("  isZ4: {0}".format(spaceicegraph.isZ4()))
    logger.debug("  defects: {0}".format(spaceicegraph.bernal_fowler_defects()))
    spaceicegraph.vector_check()
    s = "" # for yaplot

    # TSL
    # defect-defect chains
    defects = []
    for node in spaceicegraph.nodes():
        if spaceicegraph.degree(node) != 4:
            assert spaceicegraph.degree(node) < 4
            defects.append(node)
    logger.info("  Non Z4: {0}".format(defects))
    if len(defects) > 0:
        reject_count = 10
        while reject_count > 0:
            net_polar = spaceicegraph.net_polarization()
            # logger.info("  Net polarization: {0}".format(net_polar))
            if np.dot(net_polar, net_polar) < 0.05**2:
                break  # without problem
            while True:
                orig = defects[random.randint(0,len(defects)-1)]
                if spaceicegraph.out_degree(orig) == 2:
                    break
            while True:
                dest = defects[random.randint(0,len(defects)-1)]
                if spaceicegraph.in_degree(dest) == 2:
                    break
            path = shortest_path(spaceicegraph, orig, [dest,])
            if path is None:
                continue
            logger.debug("  Dipole moment = {0}".format(spaceicegraph.dipole_moment(path)))
            new_net_polar = net_polar - spaceicegraph.dipole_moment(path)*2
            if np.linalg.norm(new_net_polar) < np.linalg.norm(net_polar):
                if draw is not None:
                    s += yp.Size(0.03)
                    s += draw.draw_cell()
                    s += draw.draw_path(path)
                    s += yp.NewPage()
                spaceicegraph.invert_path(path)
                ## chk_net_polar = spaceicegraph.net_polarization()
                net_polar = new_net_polar
                logger.info("  Net polarization: {0}".format(net_polar))
                ## assert np.linalg.norm(new_net_polar - chk_net_polar) < 1e-5
                ## logger.debug("New dipole {0}".format(new_net_polar))
                ## logger.debug("Chk dipole {0}".format(chk_net_polar))
                reject_count = 10
            else:
                logger.debug("  Reject inversion")
                reject_count -= 1




    while True:
        net_polar = spaceicegraph.net_polarization()
        logger.info("  Net polarization: {0}".format(net_polar))
        if np.dot(net_polar, net_polar) < 0.2**2:
            break  # without problem
        if -1 <= net_polar[0]  <= 1 and -1 <= net_polar[1]  <= 1 and -1 <= net_polar[2]  <= 1:
            logger.info("  Gave up eliminating the polarization. (2)")
            break
        if net_polar[0] > 1.0:
            logger.debug("Depolarize +X")
            axis = np.array([+1.0, 0.0, 0.0])
        elif net_polar[0] < -1.0:
            logger.debug("Depolarize -X")
            axis = np.array([-1.0, 0.0, 0.0])
        elif net_polar[1] > 1.0:
            logger.debug("Depolarize +Y")
            axis = np.array([0.0,+1.0, 0.0])
        elif net_polar[1] < -1.0:
            logger.debug("Depolarize -Y")
            axis = np.array([0.0,-1.0, 0.0])
        elif net_polar[2] > 1.0:
            logger.debug("Depolarize +Z")
            axis = np.array([0.0, 0.0,+1.0])
        elif net_polar[2] < -1.0:
            logger.debug("Depolarize -Z")
            axis = np.array([0.0, 0.0,-1.0])
        cycle = traversing_cycle(spaceicegraph, cell, axis, draw)
        if cycle is not None:
            edges = [(cycle[i], cycle[i+1]) for i in range(len(cycle)-1)]
            if len(edges) != len(set(edges)):
                logger.debug("The cycle is entangled.")
            else:
                if draw is not None:
                    s += yp.Size(0.03)
                    s += draw.draw_cell()
                    s += draw.draw_path(cycle)
                    s += yp.NewPage()
                spaceicegraph.invert_path(cycle)
                spaceicegraph.vector_check()
    
    
    logger.debug("isZ4: {0}".format(spaceicegraph.isZ4()))
    logger.debug("defects: {0}".format(spaceicegraph.bernal_fowler_defects()))
    return s


def purge_ice_defects(icegraph):
    """
    This is faster than the method in icegraph, but
    it also polarizes the graph in the course of purging.
    """
    logger = logging.getLogger()
    while len(icegraph.bernal_fowler_defects()) > 0:
        logger.info("# of defects: {0}".format(len(icegraph.bernal_fowler_defects())))
        ins = set(icegraph.excess_in_defects())
        for out in icegraph.excess_out_defects():
            while icegraph.out_degree(out) > 2:
                path = shortest_path(icegraph, out, [ins,])
                if path is not None:
                    logger.debug("# of in defects: {0}".format(len(ins)))
                    icegraph.invert_path(path)
                    end = path[-1]
                    if icegraph.in_degree(end) == 2:
                        ins.remove(end)
                #logger.debug("IN:{0}".format(len(set(icegraph.excess_in_defects()))))
                #logger.debug("OUT:{0}".format(len(set(icegraph.excess_out_defects()))))
