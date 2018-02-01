# coding: utf-8
"""
Test suite: ring direction statistics

Hollins, G. T. Configurational statistics and the dielectric constant of ice. Proc. Phys. Soc. 84, 1001–1016 (1964).
"""

# A directed cycle is expressed as an array of True and False.
# For example,
# --> --> --> : [True, True, False]
# --> <-- --> : [Frue, False, True]
#They are isomorphic if they can be identical by rotations and inversions.

# 


import itertools as it
from countrings import countrings_nx as cr
import sys
from collections import defaultdict
import networkx as nx
from math import log
import fractions

def isomorphs(a):
    """
    return the symmetry number of cycle a
    """
    iso = set()
    aa = a + a
    for i in range(len(a)):
        part = tuple(aa[i:i+len(a)])
        iso.add(part)
    e = [not x for x in a]
    ee = e + e
    for i in range(len(a)):
        part = tuple(ee[i:i+len(a)])
        iso.add(part)
    return iso


def code(ori):
    """
    from an array of True/False to an integer.
    """
    s = 0
    for i in range(len(ori)):
        if ori[i]:
            s += 1<<i
    return s


def encode(ori):
    m = 9999
    for i in isomorphs(ori):
        c = code(i)
        if c < m:
            m = c
    return m

    
def symmetry(a):
    return len(isomorphs(a))


def freedom(a):
    """
    return possible number of molecular placements
    """
    n = 1
    for i in range(len(a)):
        if a[i-1] is a[i]:
            n *= 2
    return n


def ideal_frequency(a):
    N = len(a)
    return symmetry(a) * freedom(a)


def contains(a,b):
    """
    return True if one of the isomorphs of b is included in the set a.
    """
    return len(a & isomorphs(b)) > 0
    

def ideal_frequencies(N):
    freq = dict()
    for a in it.product((False, True), repeat=N):
        c = encode(a)
        if c not in freq:
            freq[c] = ideal_frequency(a)
    return freq


def probabilities(N):
    freq = ideal_frequencies(N)
    sum = 0
    for v in freq:
        sum += freq[v]
    prob = dict()
    for v in freq:
        prob[v] = fractions.Fraction(freq[v], sum)
    return prob


def orientations(members, digraph):
    return [digraph.has_edge(members[i-1], members[i]) for i in range(len(members))]


def hook4(lattice):
    lattice.logger.info("Hook4: Ring test.")
    graph = nx.Graph(lattice.spacegraph) #undirected
    stat = dict()
    prob = defaultdict(int)
    for n in 3,4,5,6,7,8:
        prob[n] = probabilities(n)
        stat[n] = defaultdict(int)
    for ring in cr.CountRings(graph).rings_iter(8):
        ori = orientations(ring, lattice.spacegraph)
        c   = encode(ori)
        n   = len(ring)
        stat[n][c] += 1
    #size code code(binary) Approx. Stat.
    for n in 3,4,5,6,7,8:
        fmtstr = "{{0}} {{1}} {{1:0{0}b}} {{2}} {{3:.5f}} {{4:.5f}}".format(n)
        denom = 0
        for c in stat[n]:
            denom += stat[n][c]
        if denom > 0:
            dKL = 0.0
            for c in prob[n]:
                print(fmtstr.format(n,c,prob[n][c], float(prob[n][c]), stat[n][c]/denom))
                q = stat[n][c]/denom
                p = prob[n][c]
                if q > 0.0:
                    dKL += q*(log(q) - log(p))
            dKL /= log(2.0)  #bit
            print("{1} {0} dKL[{0}-ring]".format(n,dKL))
            lattice.logger.info("  dKL[{0}-ring]={1}".format(n,dKL))
    #print(score)
    lattice.logger.info("Hook4: end.")


hooks = {4:hook4}
