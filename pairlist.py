#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import itertools
import numpy as np


def ArrangeAddress(xyz,grid,box):
    #residents in each grid cell
    residents = dict()
    for i in range(len(xyz)):
        mol = xyz[i]
        mol -= np.floor( mol / box ) * box
        address = tuple((mol / box * grid).astype(int))
        if address not in residents:
            residents[address] = set()
        residents[address].add(i)
    return residents



def pairlist(xyz,rc,box):
    #residents in each grid cell
    grid = (np.floor(box/rc)).astype(int)
    #print "START Arrange"
    residents = ArrangeAddress(xyz,grid,box)
    #print "END Arrange"

    pair = set()
    #key-value pairs in the dictionary
    donecellpair = set()
    for address in residents:
        members = residents[address]
        ix,iy,iz = address
        #neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid)%grid for j in range(-1,2)])
        for a2 in itertools.product(k[:,0],k[:,1],k[:,2]):
            if address == a2:
                #print("ISOCELL",address,npa,a2)
                for a,b in itertools.combinations(members,2):
                    pair.add((a,b))
            else:
                if a2 in residents:
                    if not (frozenset((address,a2)) in donecellpair):
                        donecellpair.add(frozenset((address,a2)))
                        #print("HETEROCELL",address,a2,donecellpair)
                        for a in members:
                            for b in residents[a2]:
                                pair.add((a,b))
    #print "PAIRLIST finished"
    return pair


#assume xyz and box are numpy.array
def pairlist_fine(xyz,rc,box):
    newpairs = []
    for i,j in pairlist(xyz,rc,box):
        moli = xyz[i]
        molj = xyz[j]
        d = moli-molj
        d -= np.floor( d/box + 0.5 )*box
        rr = np.dot(d,d)
        if rr < rc**2:
            newpairs.append((i,j,math.sqrt(rr)))
    return newpairs


def test():
    xyz = []
    for x in range(4):
        for y in range(4):
            for z in range(4):
                xyz.append(np.array((x,y,z)))
    box = np.array((4,4,4))
    rc = 1.1
    pairs = pairlist_fine(xyz,rc,box)
    for i,j,l in pairs:
        print(i,j)
        print(j,i)
    print(len(pairs))

if __name__ == "__main__":
    test()
