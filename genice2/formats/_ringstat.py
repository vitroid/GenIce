# coding: utf-8
"""
A formatter plugin to show the bond orientation statistics along the rings in an ice strucure.

Usage:
    genice2 5 -f _ringstat[max=8]

_ringstat plugin makes the statistics of bond orientations along each
HB ring and compare the distribution with that of an ideal (isolated
random) ring. The difference in the distribution is evaluated by
Kullback-Leibler divergence, d_{KL}.

A typical dKL is zero for hydrogen-disordered ices, while it is
larger than 1 for hydrogen-ordered ones like ices 2 and 9.

Ringstat analysis validates the ring-scale randomness. GenIce tool
also certifies the zero net dipole moment and Bernal-Fowler-Pauling
ice rule in terms of the validity in global and local structures.

Options:
  max=8     Size of the largest cycle to be examined.

Columns in the output:
  ring size
  code (decimal) indicating the orientations of the bonds along a cycle.
  code (binary)
  expectation (fractional) for an isolated random ring
  expectation (numerical)
  observation (fractional) in the given structure
  observation (numerical)

Finally, dKL between expectiations and observations is also shown.

% genice2 1h -r 2 2 3 -f _ringstat[max=6]

6 0 000000 64/365 0.17534 66/384 0.17188
6 1 000001 96/365 0.26301 92/384 0.23958
6 3 000011 96/365 0.26301 105/384 0.27344
6 5 000101 24/365 0.06575 21/384 0.05469
6 7 000111 48/365 0.13151 59/384 0.15365
6 9 001001 12/365 0.03288 14/384 0.03646
6 11 001011 24/365 0.06575 23/384 0.05990
6 21 010101 1/365 0.00274 4/384 0.01042
0.015526521932221845 6 dKL[6-ring]


"""

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
import numpy as np
import fractions
from math import log
import networkx as nx
from collections import defaultdict
import sys
from cycless.cycles import cycles_iter
import itertools as it

desc = {
    "ref": {
        "Hollins1964": "Hollins, G. T. Configurational statistics and the dielectric constant of ice. Proc. Phys. Soc. 84, 1001–1016 (1964)."
    },
    "brief": "Bond direction statistics.",
    "usage": __doc__,
}


# A directed cycle is expressed as an array of True and False.
# For example,
# --> --> <-- : [True, True, False]
# --> <-- --> : [True, False, True]
# They are isomorphic if they can be identical by rotations and inversions.

#


def isomorphs(a):
    """
    return the symmetry number of cycle a
    """
    iso = set()
    aa = a + a
    for i in range(len(a)):
        part = tuple(aa[i : i + len(a)])
        iso.add(part)
    e = [not x for x in a]
    ee = e + e
    for i in range(len(a)):
        part = tuple(ee[i : i + len(a)])
        iso.add(part)
    return iso


def code(ori):
    """
    from an array of True/False to an integer.
    """
    s = 0
    for i in range(len(ori)):
        if ori[i]:
            s += 1 << i
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
        if a[i - 1] is a[i]:
            n *= 2
    return n


def ideal_frequency(a):
    N = len(a)
    return symmetry(a) * freedom(a)


def contains(a, b):
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
    return [digraph.has_edge(members[i - 1], members[i]) for i in range(len(members))]


class Format(genice2.formats.Format):
    """
    _ringstat plugin makes the statistics of bond orientations along each
    HB ring and compare the distribution with that of an ideal (isolated
    random) ring. The difference in the distribution is evaluated by
    Kullback-Leibler divergence, d_{KL}.
    A typical dKL is zero for hydrogen-disordered ices, while it is
    larger than 1 for hydrogen-ordered ones like ices 2 and 9.

    Ringstat analysis validates the ring-scale randomness. GenIce tool
    also certifies the zero net dipole moment and Bernal-Fowler-Pauling
    ice rule in terms of the validity in global and local structures.

    Options:
      max=8     Size of the largest cycle to be examined.

    Columns in the output:
      ring size
      code (decimal) indicating the orientations of the bonds.
      code (binary)
      expectation (fractional) for an isolated random ring
      expectation (numerical)
      observation (fractional) in the given structure
      observation (numerical)

    dKL between expectiations and observations is also calculated.
    """

    largestring = 8

    def __init__(self, **kwargs):
        logger = getLogger()
        unknown = dict()
        for k, v in kwargs.items():
            if k == "max":
                self.largestring = int(v)
            else:
                logger.error(f"Unknown keyword: {k}")
                sys.exit(1)
        logger.info("  Largest ring: {0}.".format(self.largestring))
        super().__init__(**kwargs)

    def hooks(self):
        return {4: self.Hook4}

    @timeit
    @banner
    def Hook4(self, ice):
        "Statistics on the HBs along a ring."
        logger = getLogger()
        graph = nx.Graph(ice.digraph)  # undirected
        stat = dict()
        prob = defaultdict(int)

        # Ideal distributions
        for n in range(3, self.largestring + 1):
            prob[n] = probabilities(n)
            stat[n] = defaultdict(int)

        for ring in cycles_iter(graph, self.largestring, pos=ice.reppositions):
            ori = orientations(ring, ice.digraph)
            c = encode(ori)
            n = len(ring)
            stat[n][c] += 1

        # size code code(binary) Approx. Stat.
        logger.info(
            """
        """
        )
        s = ""
        for n in range(3, self.largestring + 1):
            fmtstr = (
                "{{0}} {{1}} {{1:0{0}b}} {{2}} {{3:.5f}} {{5}}/{{6}} {{4:.5f}} ".format(
                    n
                )
            )
            denom = 0
            for c in stat[n]:
                denom += stat[n][c]
            if denom > 0:
                dKL = 0.0
                for c in prob[n]:
                    s += (
                        fmtstr.format(
                            n,
                            c,
                            prob[n][c],
                            float(prob[n][c]),
                            stat[n][c] / denom,
                            stat[n][c],
                            denom,
                        )
                        + "\n"
                    )
                    q = stat[n][c] / denom
                    p = prob[n][c]
                    if q > 0.0:
                        dKL += q * (log(q) - log(p))
                dKL /= log(2.0)  # bit
                s += "{1} {0} dKL[{0}-ring]\n".format(n, dKL)
                logger.info("  dKL[{0}-ring]={1}".format(n, dKL))
        self.output = s
