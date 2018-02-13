#!/usr/bin/env python3

import sys
import numpy as np
import itertools as it

def lookup_and_replace(v1,pvecs):
    L1 = sum([x**2 for x in v1])
    for i,v2 in enumerate(pvecs):
        L2 = sum([x**2 for x in v2])
        ip = np.dot(v1,v2)
        if ip**2 == L1*L2:
            if L1 < L2:
                pvecs[i] = v1
            return
    pvecs.append(v1)
    return


def primevectors():
    pvecs = []
    for v1 in it.product(range(10,-10-1,-1),repeat=3):
        print(v1)
        if v1 == (0,0,0):
            continue
        lookup_and_replace(np.array(v1).astype(float),pvecs)
    return pvecs
                
            
tolerance = 0.02
cell = np.array([float(x) for x in sys.argv[1:4]])

#Reference vector
x = np.array([1,1,1])
Lx = np.dot(x*cell,x*cell)
pvecs = np.array(primevectors())
#print(pvecs)
#inner product of prime vectors and x
ipx = np.dot(pvecs[:]*cell,x*cell)
#print(ipx)
#length(pvecs)**2
ipp = np.sum(pvecs*pvecs*cell*cell, axis=1)
#ipp[:] = np.linalg.norm(pvecs[:,])
#print(ipp)
#vectors orthogonal to x
cond = ipx*ipx < ipp*Lx * tolerance**2
ys = pvecs[cond]
ipp = ipp[cond]
for i in range(ys.shape[0]):
    vi = ys[i]
    for j in range(i+1,ys.shape[0]):
        vj = ys[j]
        vij = np.dot(vi*cell,vj*cell)
        if vij*vij < ipp[i]*ipp[j]*tolerance**2:
            v0 = (ipp[i]*ipp[j]*Lx)**0.5
            v1 = np.dot(np.cross(x*cell,vi*cell),vj*cell)
            if v1 < 0:
                v1 = -v1
                vj = -vj
            print(v1,v1/v0,end="")
            s = " {0} {1} {2} ".format(*x)
            s += "{0:.0f} {1:.0f} {2:.0f} ".format(*vi)
            s += "{0:.0f} {1:.0f} {2:.0f} ".format(*vj)
            print(s)

    
