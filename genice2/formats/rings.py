# coding: utf-8
"""
Show rings in Yaplot format, defined in https://github.com/vitroid/Yaplot

Usage:
    genice2 III -f rings > 3.rings.yap
    genice2 III -f rings[5] > output.yap    # up to 5-membered rings.
"""

desc = { "ref": { "MBO2007": "Matsumoto, M., Baba, A. & Ohmine, I. Topological building blocks of hydrogen bond network in water. J. Chem. Phys. 127, 134504 (2007).",
                  "CountRings": "https://github.com/vitroid/CountRings",
                  "Yaplot":     "https://github.com/vitroid/Yaplot" },
         "brief": "Show rings in Yaplot.",
         "usage": __doc__
         }



import sys
from collections import defaultdict
import numpy as np
import networkx as nx
import yaplotlib as yp
from logging import getLogger

from genice2 import rigid
import genice2.formats
from countrings import countrings_nx as cr



def face(center, rpos):
    pos = rpos + center
    n = rpos.shape[0]
    s = yp.Color(n)
    s += yp.Layer(n)
    s += yp.Polygon(pos)
    return s


class Format(genice2.formats.Format):
    largestring=8

    def __init__(self, **kwargs):
        logger = getLogger()
        unknown = dict()
        for k, v in kwargs.items():
            try:
                self.largestring=int(k)
            except:
                logger.error("Argument must be a positive integer.")
                sys.exit(1)
        logger.info("  Largest ring: {0}.".format(self.largestring))
        super().__init__()


    def hooks(self):
        return {2:self.hook2}


    def hook2(self, ice):
        logger = getLogger()
        logger.info("Hook2: Show rings in Yaplot format.")
        # copied from svg_poly
        graph = nx.Graph(ice.graph) #undirected
        cellmat = ice.repcell.mat
        s = ""
        s += yp.Layer(2)
        s += yp.Color(0)
        for i,j in graph.edges():
            pi, pj = ice.reppositions[i], ice.reppositions[j]
            d = pj - pi
            d -= np.floor(d+0.5)
            s += yp.Line(pi @ cellmat, (pi+d) @ cellmat)
        for ring in cr.CountRings(graph, pos=ice.reppositions).rings_iter(self.largestring):
            deltas = np.zeros((len(ring),3))
            d2 = np.zeros(3)
            for k,i in enumerate(ring):
                d = ice.reppositions[i] - ice.reppositions[ring[0]]
                d -= np.floor(d+0.5)
                deltas[k] = d
            comofs = np.sum(deltas, axis=0) / len(ring)
            deltas -= comofs
            com = ice.reppositions[ring[0]] + comofs
            com -= np.floor(com)
            # rel to abs
            com    = np.dot(com,    cellmat)
            deltas = np.dot(deltas, cellmat)
            s += face(com,deltas)
        self.output = s + yp.NewPage()
        logger.info("Hook2: end.")
