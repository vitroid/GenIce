# coding: utf-8
"""
Test suite: Kirkwood G factor
"""

import numpy as np

#It should be expressed as a function of distance. 

def hook5(lattice):
    lattice.logger.info("Hook5: Kirkwood G.")
    zvec = np.array([0.,0.,1.])
    dipoles = np.dot(zvec, lattice.rotmatrices)
    n = dipoles.shape[0]
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
    lattice.logger.info("Hook5: end.")


hooks = {5:hook5}
