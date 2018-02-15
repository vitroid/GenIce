# coding: utf-8
"""
Yaplot format.
defined in https://github.com/vitroid/Yaplot
"""

from collections import defaultdict
import numpy as np
import networkx as nx

from genice import rigid
from genice import yaplotlib as yp
from countrings import countrings_nx as cr




def face(center, rpos):
    pos = rpos + center
    n = rpos.shape[0]
    s = yp.Color(n)
    s += yp.Layer(n)
    s += yp.Polygon(pos)
    return s



def hook4(lattice):
    lattice.logger.info("Hook4: A. Output depolarization process in Yaplot format.")
    print(lattice.yapresult, end="")
    lattice.logger.info("Hook4: B. Polyhedral expression.")
    # copied from svg_poly
    graph = nx.Graph(lattice.spacegraph) #undirected
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
    lattice.logger.info("Hook4: end.")


def hook6(lattice):
    global nwateratoms
    lattice.logger.info("Hook6: Output water molecules in Yaplot format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    # prepare the reverse dict
    waters = defaultdict(dict)
    for atom in lattice.atoms:
        resno, resname, atomname, position, order = atom
        if "O" in atomname:
            waters[order]["O"] = position
        elif "H" in atomname:
            if "H0" not in waters[order]:
                waters[order]["H0"] = position
            else:
                waters[order]["H1"] = position
    s = ""
    s += yp.Color(3)
    for order, water in waters.items():
        O = water["O"]        
        H0 = water["H0"]        
        H1 = water["H1"]        
        s += yp.Layer(4)
        s += yp.Color(3)
        s += yp.Size(0.03)
        s += yp.Circle(O)
        s += yp.Line(O,H0)
        s += yp.Line(O,H1)
        s += yp.Size(0.01)
        s += yp.Circle(H0)
        s += yp.Circle(H1)
        s += yp.Color(2)
        s += yp.Layer(1)
        s += yp.Text(O, "{0}".format(order))
    s += yp.Layer(3)
    s += yp.Color(4)
    s += yp.ArrowType(1)
    s += yp.Size(0.03)
    for i,j in lattice.spacegraph.edges(data=False):
        if i in waters and j in waters:  # edge may connect to the dopant
            O = waters[j]["O"]
            H0 = waters[i]["H0"]
            H1 = waters[i]["H1"]
            d0 = H0 - O
            d1 = H1 - O
            rr0 = np.dot(d0,d0)
            rr1 = np.dot(d1,d1)
            if rr0 < rr1 and rr0 < 0.27**2:
                s += yp.Arrow(H0,O)
            if rr1 < rr0 and rr1 < 0.245**2:
                s += yp.Arrow(H1,O)
    print(s, end="")
    nwateratoms = len(lattice.atoms)
    lattice.logger.info("Hook6: end.")


def hook7(lattice):
    global nwateratoms
    lattice.logger.info("Hook7: Output water molecules in Yaplot format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    gatoms = lattice.atoms[nwateratoms:]
    palettes = dict()
    s = ""
    s += yp.Layer(4)
    s += yp.ArrowType(1)
    H = []
    O  = ""
    for atom in gatoms:
        resno, resname, atomname, position, order = atom
        if atomname in palettes:
            pal = palettes[atomname]
        else:
            pal = 4 + len(palettes)
            palettes[atomname] = pal
        s += yp.Color(pal)
        s += yp.Size(0.04)
        s += yp.Circle(position)
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s)
    lattice.logger.info("Hook7: end.")
    

hooks = {7:hook7, 6:hook6, 4:hook4}
