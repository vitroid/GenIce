#!/usr/bin/env python
# -*- coding: utf-8 -*-

####Note: xyz are in coord relative to the cell.


from __future__ import print_function
import math
import itertools
import numpy as np
import logging

def ArrangeAddress(xyz,grid):
    #residents in each grid cell
    residents = dict()
    for i in range(len(xyz)):
        mol = xyz[i]
        mol -= np.floor( mol )
        address = tuple((mol * grid).astype(int))
        if address not in residents:
            residents[address] = set()
        residents[address].add(i)
    return residents



def pairlist(xyz,grid):
    #print "START Arrange"
    residents = ArrangeAddress(xyz,grid)
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


def pairlist_hetero(xyz,xyz2,grid):
    #print "START Arrange"
    residents  = ArrangeAddress(xyz,grid)
    residents2 = ArrangeAddress(xyz2,grid)
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
            if a2 in residents2:
                if not ((address,a2) in donecellpair):
                        donecellpair.add((address,a2))
                        #print("HETEROCELL",address,a2,donecellpair)
                        for a in members:
                            for b in residents2[a2]:
                                pair.add((a,b))
    #print "PAIRLIST finished"
    return pair


#assume xyz and box are numpy.array
def pairlist_fine(xyz,rc,cell,grid,distance=True):
    newpairs = []
    for i,j in pairlist(xyz,grid):
        moli = xyz[i]
        molj = xyz[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
            
        if rr < rc**2:
            if distance:
                newpairs.append((i,j,math.sqrt(rr)))
            else:
                newpairs.append((i,j))
    return newpairs


#assume xyz and box are numpy.array
def pairlist_fine_hetero(xyz,xyz2,rc,cell,grid,distance=True):
    newpairs = []
    for i,j in pairlist_hetero(xyz,xyz2,grid):
        moli = xyz[i]
        molj = xyz2[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
            
        if rr < rc**2:
            if distance:
                newpairs.append((i,j,math.sqrt(rr)))
            else:
                newpairs.append((i,j))
    return newpairs


#
def determine_grid(cell, radius):
    logger = logging.getLogger()
    ct = cell.transpose()
    a = ct[0]
    b = ct[1]
    c = ct[2]
    al = np.linalg.norm(a)   #vector length
    bl = np.linalg.norm(b)
    cl = np.linalg.norm(c)
    ae = a / al              #unit vectors
    be = b / bl
    ce = c / cl
    ad = np.dot(ae,np.cross(be,ce)) #distance to the bc plane
    bd = np.dot(be,np.cross(ce,ae))
    cd = np.dot(ce,np.cross(ae,be))
    ax = radius / ad        # required length of a vector to contain a sphere of radius 
    bx = radius / bd
    cx = radius / cd
    gf = np.array([al/ax, bl/bx, cl/cx])  # required number of grid cells
    #print(cell,radius,gf)
    #import sys
    #sys.exit(1)
    return np.floor(gf).astype(int)



def test():
    xyz = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz.append(np.array((x/100.+1,y/100.+1,z/100.+1)) / 4.)
    xyz2 = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz2.append(np.array((x/100.+2,y/100.+1,z/100.+1)) / 4.)
    box = np.diag((4,4,4))
    rc = 1.000000001
    grid = determine_grid(box,rc)
    pairs = pairlist_fine_hetero(xyz,xyz2,rc,box,grid)
    for i,j,l in pairs:
        print(i,j,l)
    print(len(pairs))

if __name__ == "__main__":
    test()
