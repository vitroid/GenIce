# coding: utf-8
"""
"""

import numpy as np
# import itertools as it
from countrings import countrings_nx as cr
# import sys
# from collections import defaultdict
import networkx as nx
# from math import log
# import fractions
import re


offset = np.zeros(3)


# copied from povray.py

def Block(name, content):
    return " " + name + " { " + content + " } "

def Vector(v):
    return "<{0:.3f}, {1:.3f}, {2:.3f}> ".format(*v)

def Juxtapose(v):
    return ",".join(v)

def Atom(atomtype, pos):
    return Block( "sphere", Juxtapose( [Vector(pos), "R{0}".format(atomtype) ] ) + Block( "material", "MAT{0}".format(atomtype) ) ) + "\n"
    
def Bond(bondtype, pos1, pos2):
    return Block( "cylinder", Juxtapose( [Vector(pos1), Vector(pos2), "R{0}".format(bondtype)] ) + Block( "material", "MAT{0}".format(bondtype) ) ) + "\n"

# end copy


def SmoothTriangle(facetype, pos):
    vecs = [Vector(p) for p in pos]
    return Block( "smooth_triangle", Juxtapose( vecs ) + Block( "material", "MAT{0}".format(facetype) ) ) + "\n"



def face(center, rpos):
    n = rpos.shape[0]
    s = "// {0}\n".format(n)
    # draw skeleton
    for i in range(n):
        p1 = center + rpos[i-1]
        p2 = center + rpos[i]
        s += Bond("HB", p1, p2)
    # normalize relative vectors
    normalized = np.zeros_like(rpos)
    for i in range(n):
        normalized[i] = rpos[i] / np.linalg.norm(rpos[i])
    #normal for each triangle
    normals = np.zeros_like(rpos)
    for i in range(n):
        normals[i] = np.cross(normalized[i-1], normalized[i])
        normals[i] /= np.linalg.norm(normals[i])
    # central normal
    c_normal = np.sum(normals, axis=0) / n
    # normal for each vertex
    v_normal = np.zeros_like(rpos)
    for i in range(n):
        v_normal[i-1] = normals[i-1] + normals[i]
        v_normal[i-1] /= np.linalg.norm(v_normal[i-1])
    for i in range(n):
        s += SmoothTriangle("r{0}".format(n), [center, c_normal, center+rpos[i-1], v_normal[i-1], center+rpos[i], v_normal[i] ])
    return s


def hook4(lattice):
    global offset
    lattice.logger.info("Hook4: Povray (polyhedral expressions).")
    graph = nx.Graph(lattice.spacegraph) #undirected
    cellmat = lattice.repcell.mat
    queue = []
    for ring in cr.CountRings(graph).rings_iter(8):
        deltas = np.zeros((len(ring),3))
        for k,i in enumerate(ring):
            d = lattice.reppositions[i] - lattice.reppositions[ring[0]]
            d -= np.floor(d+0.5)
            deltas[k] = d
        comofs = np.sum(deltas, axis=0) / len(ring)
        deltas -= comofs
        com = lattice.reppositions[ring[0]] + comofs + offset
        com -= np.floor(com)
        # rel to abs
        com    = np.dot(com,    cellmat)
        deltas = np.dot(deltas, cellmat)
        queue.append((com, deltas))
    s = '#include "default.inc"' + "\n"
    s += "union {\n"
    # s = '//' + "\n//".join(lattice.doc) + "\n" + s
    for com, deltas in queue:
        s += face(com, deltas)
    s += "  translate " + Vector( -(cellmat[0,:]+cellmat[1,:]+cellmat[2,:])/2 ) + "\n}\n"
    print(s)
    lattice.logger.info("Hook4: end.")


def argparser(arg):
    global offset
    assert re.match("^[-+0-9,.]+$", arg) is not None, "Argument must be three floating points separated by commas."
    offset = np.array([float(x) for x in arg.split(",")]).reshape(3)
        

hooks = {4:hook4}
