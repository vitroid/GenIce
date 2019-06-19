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
import numpy as np

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


def spanning(ring, rpos):
    dsum = np.zeros(3)

    for k in range(len(ring)):
        d = rpos[ring[k-1]] - rpos[ring[k]]
        d -= np.floor(d+0.5)
        dsum += d

    return not np.all(np.absolute(dsum) < 1e-5)
    



def hook4(lattice):
    lattice.logger.info("Hook4: Statistics on the HBs along a ring.")
    graph = nx.Graph(lattice.spacegraph) #undirected
    stat = dict()
    prob = defaultdict(int)

    # Ideal distributions
    for n in range(3, lattice.largestring+1):
        prob[n] = probabilities(n)
        stat[n] = defaultdict(int)

    for ring in cr.CountRings(graph).rings_iter(lattice.largestring):
        if not spanning(ring, lattice.reppositions):
            ori = orientations(ring, lattice.spacegraph)
            c   = encode(ori)
            n   = len(ring)
            stat[n][c] += 1

    #size code code(binary) Approx. Stat.
    lattice.logger.info("""
    _ringstat plugin makes the statistics of bond orientations along each
    HB ring and compare the distribution with that of an ideal (isolated
    random) ring. The difference in the distribution is evaluated by
    Kullback-Leibler # divergence, d_{KL}.
    A typical dKL is zero for hydrogen-disordered ices, while it is
    larger than 1 for hydrogen-ordered ones like ices 2 and 9.

    Ringstat analysis validates the ring-scale randomness. GenIce tool
    also certifies the zero net dipole moment and Bernal-Fowler-Pauling
    ice rule in terms of the validity in global and local structures.

    Columns in the output:
      ring size
      code (decimal) indicating the orientations of the bonds.
      code (binary)
      expectation (fractional) for an isolated random ring
      expectation (numerical)
      observation (fractional) in the given structure
      observation (numerical)

    dKL between expectiations and observations is also calculated.
    """)
    for n in range(3, lattice.largestring+1):
        fmtstr = "{{0}} {{1}} {{1:0{0}b}} {{2}} {{3:.5f}} {{5}}/{{6}} {{4:.5f}} ".format(n)
        denom = 0
        for c in stat[n]:
            denom += stat[n][c]
        if denom > 0:
            dKL = 0.0
            for c in prob[n]:
                print(fmtstr.format(n,c,prob[n][c], float(prob[n][c]), stat[n][c]/denom, stat[n][c], denom))
                q = stat[n][c]/denom
                p = prob[n][c]
                if q > 0.0:
                    dKL += q*(log(q) - log(p))
            dKL /= log(2.0)  #bit
            print("{1} {0} dKL[{0}-ring]".format(n,dKL))
            lattice.logger.info("  dKL[{0}-ring]={1}".format(n,dKL))
    #print(score)
    lattice.logger.info("Hook4: end.")


# argparser
def hook0(lattice, arg):
    lattice.logger.info("Hook0: ArgParser.")

    if arg == "":
        lattice.largestring=8
    else:
        try:
            lattice.largestring=int(arg)
        except:
            lattice.logger.error("Argument must be a positive integer.")
            sys.exit(1)

    lattice.logger.info("  Largest ring: {0}.".format(lattice.largestring))
    lattice.logger.info("Hook0: end.")


hooks = {4:hook4, 0:hook0}
