# coding: utf-8
"""
Show rings in Yaplot format.
defined in https://github.com/vitroid/Yaplot
"""

from collections import defaultdict
import numpy as np
import networkx as nx
import yaplotlib as yp

from genice import rigid
from countrings import countrings_nx as cr




def face(center, rpos):
    pos = rpos + center
    n = rpos.shape[0]
    s = yp.Color(n)
    s += yp.Layer(n)
    s += yp.Polygon(pos)
    return s



def hook2(lattice):
    lattice.logger.info("Hook2: Show rings in Yaplot format.")
    # copied from svg_poly
    graph = nx.Graph(lattice.graph) #undirected
    cellmat = lattice.repcell.mat
    for ring in cr.CountRings(graph).rings_iter(8):
        deltas = np.zeros((len(ring),3))
        d2 = np.zeros(3)
        for k,i in enumerate(ring):
            d = lattice.reppositions[i] - lattice.reppositions[ring[0]]
            d -= np.floor(d+0.5)
            deltas[k] = d
            dd = lattice.reppositions[ring[k]] - lattice.reppositions[ring[k-1]]
            dd -= np.floor(dd+0.5)
            d2 += dd
        # d2 must be zeros
        if np.all(np.absolute(d2) < 1e-5):
            comofs = np.sum(deltas, axis=0) / len(ring)
            deltas -= comofs
            com = lattice.reppositions[ring[0]] + comofs
            com -= np.floor(com)
            # rel to abs
            com    = np.dot(com,    cellmat)
            deltas = np.dot(deltas, cellmat)
            print(face(com,deltas), end="")
    print()
    lattice.logger.info("Hook2: end.")


    

hooks = {2:hook2}
