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


def estimate_density(atoms, cell, bondlen):
    logger = logging.getLogger()
    dmin = shortest_distance(atoms, cell)
    logger.debug(dmin)
    scale = bondlen / dmin
    return 18/6.022e23*len(atoms) / (np.linalg.det(cell)*1e-24 * scale**3)

    


def is_zero(v):
    return np.dot(v,v) < 1e-10


def equivalents(v):
    if is_zero(v[0]-0.5):
        v[0] -= 1.0
    if is_zero(v[1]-0.5):
        v[1] -= 1.0
    if is_zero(v[2]-0.5):
        v[2] -= 1.0
    img = [[0.],[0.],[0.]]
    if is_zero(v[0] + 0.5):
        img[0] = [0., 1.]
    if is_zero(v[1] + 0.5):
        img[1] = [0., 1.]
    if is_zero(v[2] + 0.5):
        img[2] = [0., 1.]
    return v, img


def tetrahedra(pairs, rc, coord, cell):
    logger = logging.getLogger()
    vertices = list(set([v for v in lg.flatten(pairs)]))
    logger.debug(vertices)
    neiv = dict()
    neid = dict()
    for v in vertices:
        neiv[v] = []
        neid[v] = []
    for i,j in pairs:
        d = coord[j] - coord[i]
        d -= np.floor(d + 0.5)
        d, img = equivalents(d)
        logger.debug((d,img))
        for x in img[0]:
            for y in img[1]:
                for z in img[2]:
                    dd = d + np.array((x,y,z))
                    neid[i].append(dd)
                    neiv[i].append(j)
                    neid[j].append(-dd)
                    neiv[j].append(i)
    for v in vertices:
        logger.debug(len(neiv[v]))
        for i,j,k in it.combinations(range(len(neiv[v])), 3):
            vi, vj, vk = neiv[v][i], neiv[v][j], neiv[v][k]
            if vi < v or vj < v or vk < v:
                continue
            di, dj, dk = neid[v][i], neid[v][j], neid[v][k]
            dij = np.dot(di - dj, cell)
            djk = np.dot(dj - dk, cell)
            dki = np.dot(dk - di, cell)
            if np.dot(dij,dij) < rc**2:
                if np.dot(djk,djk) < rc**2:
                    if np.dot(dki,dki) < rc**2:
                        logger.debug((v,vi,vj,vk))
                        yield (v,vi,vj,vk), (coord[v], di,dj,dk)


def toWater(coord, cell):
    logger = logging.getLogger()
    dmin  = shortest_distance(coord,cell)
    pairs = [v for v in pl.pairlist_crude(coord, dmin*1.4, cell, distance=False)]
    for vtet, dtet in tetrahedra(pairs, dmin*1.4, coord, cell):
        p = dtet[0] + (dtet[1] + dtet[2] + dtet[3])/4
        p -= np.floor( p )
        yield p

    
    
