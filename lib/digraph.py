#!/usr/bin/env python3
#coding: utf-8
#input: coordinate of the nodes, the digraph obeying the ice rule.
#output: the digraph with zero net dipole.

import sys
import math
import networkx
import random
import numpy as np
import logging

class IceGraph(networkx.DiGraph):
    #def __init__(data=None):
    #    super(IceGraph, self).__init__(data)

    def register_pairs(self,pairs):
        self.clear()
        for pair in pairs:
            x,y = pair[0:2]
            self.add_edge(x,y)

            
    def _goahead(self,node,marked,order):
        while not marked[node]:
            marked[node] = True
            order.append(node)
            nei = self.neighbors(node)
            next = random.randint(0,1)
            if len(nei) != 2:
                logging.getLogger().error("Dangling bond: {0} {1}".format(node,nei))
            node = nei[next]
        #a cyclic path is found.
        #trim the uncyclic part
        while order[0] != node:
            order.pop(0)
        #add the node at the last again
        order.append(node)
        return order


    def homodromiccycle(self):
        """
        Randomly select a homodromic cycle
        """
        marked = [False] * self.number_of_nodes()
        order = []
        node = random.randint(0,self.number_of_nodes()-1)
        return self._goahead(node,marked,order)


    def isZ4(self):
        """
        Reply whether all the vertices have four neighbors or not.
        """
        good = True
        undir = self.to_undirected()
        for node in range(undir.number_of_nodes()):
            if len(undir.neighbors(node)) != 4:
                good = False
        return good


    def purgedefects(self, defects):
        d = defects[0]
        if self.in_degree(d) == 2 and self.out_degree(d) == 2:
            defects.pop(0)
            return
        if self.in_degree(d) > 2:
            nodes = self.predecessors(d)
            i = random.randint(0,len(nodes)-1)
            node = nodes[i]
            self.remove_edge(node,d)
            self.add_edge(d,node)
            defects.append(node)
        if self.out_degree(d) > 2:
            nodes = self.successors(d)
            i = random.randint(0,len(nodes)-1)
            node = nodes[i]
            self.remove_edge(d,node)
            self.add_edge(node,d)
            defects.append(node)

            
    def defects(self):
        """
        Reply the list of defective vertices.
        """
        defects = []
        for i in range(self.number_of_nodes()):
            if self.in_degree(i) != 2 or self.out_degree(i) != 2:
                defects.append(i)
            if self.degree(i) != 4:
                logger = logging.getLogger()
                logger.error("Non Z4 vertex: {0} {1} {2} {3}".format(i,self.degree(i),self.successors(i),self.predecessors(i)))
        return defects

    def purge_ice_defects(self):
        logger = logging.getLogger()
        if not self.isZ4():
            logger.error("Some water molecules do not have four HBs.")
            sys.exit(1)
        defects = self.defects()
        while len(defects)>0:
            self.purgedefects(defects)
        if len(self.defects()) != 0:
            logger.error("Some water molecules do not obey the ice rule.")
            sys.exit(1)


class SpaceIceGraph(IceGraph):
    """
    Digraph with geometrical info
    """
    def __init__(self, data=None, coord=None, pbc=True):
        super(SpaceIceGraph, self).__init__(data)
        if coord is not None:
            self.add_vectors(coord, pbc)
            
    def add_vectors(self, coord, pbc=True):
        """
        add vector attributes to each edge
        """
        for i,j,k in self.edges_iter(data=True):
            vec = coord[j] - coord[i]
            if pbc:
                vec -= np.floor(vec + 0.5)
            k["vector"] = vec  #each edge has "vector" attribute
        
    def dipole_of_a_cycle(self, order):
        """
        normally zero.
        Non-zero when it goes across the cell.
        """
        delta = np.zeros(3)
        for i in range(len(order)-1):
            delta += self.get_edge_data(order[i],order[i+1])["vector"]
        return delta
            
    def invert_edge(self,from_,to_):
        """
        also invert the attribute vector
        """
        v = self.get_edge_data(from_,to_)["vector"]
        self.remove_edge(from_,to_)
        self.add_edge(to_,from_,vector=-v)


    def net_polarization(self):
        dipole = np.zeros(3)
        for i,j,k in self.edges_iter(data=True):
            dipole += k["vector"]
        return dipole


    #This is too slow for a big system.  Improve it.
    def depolarize(self):
        """
        It assumes vector attribute is set
        """
        logger = logging.getLogger()
        dipole = self.net_polarization()
        logger.debug("Initial dipole: {0}".format(dipole))
        s0 = np.dot(dipole, dipole)
        #In the following calculations, there must be error allowances.
        while s0 > 0.1:
            path = self.homodromiccycle()
            pathdipole = self.dipole_of_a_cycle(path)
            newdipole = dipole - 2.0 * pathdipole
            logger.debug("Updated dipole: {0}".format(newdipole))
            #logger.debug("Debugged dipole: {0}".format(self.polarization()))
            s1 = np.dot(newdipole, newdipole)
            if s1 < s0:
                #accept the inversion
                for i in range(len(path)-1):
                    f = path[i]
                    t = path[i+1]
                    self.invert_edge(f,t) #also invert the attribute vector
                s0 = s1
                dipole -= 2.0 * pathdipole
                logger.debug("Score: {0}".format(s0))
