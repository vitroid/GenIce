# coding: utf-8
"""
Hydrogen bond network in @NGPH format.
"""

import numpy as np


def kirkwood_g(lattice):
    n = len(lattice.spacegraph.coord)
    v1 = np.array([lattice.spacegraph[node][lattice.spacegraph.neighbors(node)[0]]['vector'] for node in range(n)])
    v2 = np.array([lattice.spacegraph[node][lattice.spacegraph.neighbors(node)[1]]['vector'] for node in range(n)])
    z = v1+v2
    z = np.dot(z, lattice.repcell.mat) # fractional to absolute coord
    zL = np.linalg.norm(z, axis=1)  # normalize
    for i in range(n):
        z[i] /= zL[i]
    dipoles = z
    lattice.logger.info("  Inner products of the dipoles.")
    # inner products of dipole pairs
    ip = np.zeros((n,n,3))
    for i in (0,1,2): # x,y,z
        x = dipoles[:,i]
        xT = x.reshape((n,1))
        ip[:,:,i] = x*xT  #broadcast product
    ip = np.sum(ip, axis=2).reshape(n*n)
    lattice.logger.info("  Pair distances.")
    # pair distances
    coord = lattice.reppositions
    delta = np.zeros((n,n,3))
    for i in (0,1,2):
        x = coord[:,i]
        xT = x.reshape((n,1))
        d = x - xT
        d -= np.floor( d + 0.5 )  #wrap
        delta[:,:,i] = d
    delta = np.dot(delta, lattice.repcell.mat)
    delta = delta*delta
    delta = np.sum(delta, axis=2)
    delta = np.sqrt(delta).reshape(n*n)
    # print(ip)
    # print(delta)
    # print(lattice.repcell.mat)
    lattice.logger.info("  G(r)")
    mx = np.max(delta)
    x = 0.04
    while x < mx:
        filter = ip[delta<x]
        print("{0:.5f} {1:.5f}".format(x, np.sum(filter)/n))
        x += 0.08
    



def hook4(lattice):
    lattice.logger.info("Hook4: Invert homodromic rings.")
    # unfix all edges.
    for i,j in lattice.spacegraph.edges(data=False):
        lattice.spacegraph[i][j]["fixed"] = False
    
    mon = 1
    for i in range(len(lattice.reppositions)):
        cyc = lattice.spacegraph.homodromiccycle()
        lattice.spacegraph.invert_path(cyc)
        if i == mon:
            print("# {0}".format(mon))
            kirkwood_g(lattice)
            mon += mon
            print("")
    lattice.logger.info("Hook4: end.")

hooks = {4:hook4}
