#!/usr/bin/env python3

#from Frank-Kasper (Tetrahedrally close-packed) alloy to its dual.

#Standard libs
import itertools as it
import logging
from collections import defaultdict

#non-standard libs
import numpy as np

#genice libs
from genice import pairlist as pl
from genice import libgenice as lg


def shortest_distance(atoms, cell):
    dmin = 1e99
    for a1,a2 in it.combinations(atoms,2):
        d = a1-a2
        d -= np.floor(d + 0.5)
        dv = np.dot(d, cell)
        dd = np.dot(dv,dv)
        if dd < dmin:
            dmin = dd
    return dmin**0.5


def is_zero(v):
    return np.dot(v,v) < 1e-10

def tetrahedra(pairs, coord, cell):
    logger = logging.getLogger()
    vertices = set(lg.flatten(pairs))
    logger.debug(vertices)
    nei = dict()
    for v in vertices:
        nei[v] = dict()
    for i,j in pairs:
        d = coord[j] - coord[i]
        d -= np.floor( d + 0.5 )
        nei[i][j] = d
        nei[j][i] = -d
    for v in vertices:
        logger.debug((v,len(nei[v])))
    for v in vertices:
        neis = []
        for i in nei[v]:
            if v<i:
                neis.append(i)
        logger.debug((v,neis))
        for i,j,k in it.combinations(neis, 3):
            if j in nei[i]:
                if k in nei[j]:
                    if i in nei[k]:
                        if ( is_zero(nei[v][i]+nei[i][j]+nei[j][v]) and
                             is_zero(nei[v][j]+nei[j][k]+nei[k][v]) and
                             is_zero(nei[v][k]+nei[k][i]+nei[i][v]) ):
                            logger.debug((i,j,k))
                
        


def toWater(coord, cell):
    logger = logging.getLogger()
    logger.debug(coord)
    dmin  = shortest_distance(coord,cell)
    pairs = [v for v in pl.pairlist_crude(coord, dmin*1.4, cell, distance=False)]
    tet   = tetrahedra(pairs, coord, cell)
    
    
    
