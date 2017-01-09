#!/usr/bin/env python3
#coding: utf-8
#input: coordinate of the nodes, the digraph obeying the ice rule.
#output: the digraph with zero net dipole.

import sys
import math
import networkx
import random
import numpy
import logging

class MyDiGraph(networkx.DiGraph):
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
        delta = numpy.zeros(3)
        for i in range(len(order)-1):
            delta += self.get_edge_data(order[i],order[i+1])["vector"]
        return order, delta


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
