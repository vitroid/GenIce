#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Making a dual from a Frank-Kasper (Tetrahedrally close-packed) alloy structure.
"""

# Standard libs
import itertools as it
import logging
from collections import defaultdict

# non-standard libs
import numpy as np
import pairlist as pl


def shortest_distance(atoms, cell):
    logger = logging.getLogger()
    dmin = 1e99
    for a1, a2 in it.combinations(atoms, 2):
        d = a1 - a2
        d -= np.floor(d + 0.5)
        dv = d @ cell
        dd = dv @ dv
        if dd < dmin:
            dmin = dd
    logger.debug("shortest_distance: {0}".format(dmin**0.5))
    return dmin**0.5


def estimate_density(atoms, cell, bondlen):
    dmin = shortest_distance(atoms, cell)
    scale = bondlen / dmin
    return 18 / 6.022e23 * len(atoms) / \
        (np.linalg.det(cell) * 1e-24 * scale**3)


def is_zero(v):
    return v @ v < 1e-10


def equivalents(v, cell, rc):
    """
    yield a set of vectors pointing to the image of the original point v.
    """
    origin = v.copy()
    img = [[0., 1.], [0., 1.], [0., 1.]]
    for d in range(3):
        if origin[d] > 0.0:
            origin[d] -= 1.0
    for x in img[0]:
        for y in img[1]:
            for z in img[2]:
                d = origin + np.array([x, y, z])
                r = d @ cell
                if r @ r < rc**2:
                    yield d


def adjacency_vectors(pairs, rc, coord, cell):
    logger = logging.getLogger()
    vertices = list(set([v for pair in pairs for v in pair]))
    adjv = dict()
    adjd = dict()
    for v in vertices:
        adjv[v] = []
        adjd[v] = []
    for i, j in pairs:
        d = coord[j] - coord[i]
        d -= np.floor(d + 0.5)
        for dd in equivalents(d, cell, rc):
            adjd[i].append(dd)
            adjv[i].append(j)
            adjd[j].append(-dd)
            adjv[j].append(i)
    return vertices, adjv, adjd


def tetrahedra(pairs, rc, coord, cell):
    logger = logging.getLogger()
    vertices, adjv, adjd = adjacency_vectors(pairs, rc, coord, cell)
    for v in vertices:
        logger.debug(len(adjv[v]))
        for i, j, k in it.combinations(range(len(adjv[v])), 3):
            vi, vj, vk = adjv[v][i], adjv[v][j], adjv[v][k]
            if vi < v or vj < v or vk < v:
                continue
            di, dj, dk = adjd[v][i], adjd[v][j], adjd[v][k]
            dij = (di - dj) @ cell
            djk = (dj - dk) @ cell
            dki = (dk - di) @ cell
            if dij @ dij < rc**2:
                if djk @ djk < rc**2:
                    if dki @ dki < rc**2:
                        logger.debug((v, vi, vj, vk))
                        yield (v, vi, vj, vk), (coord[v], di, dj, dk)


def toWater(coord, cell, tolerance=1.4):
    """
    cell: cell matrix (first row == a)
    coord: relative coord of the cage position
    """
    logger = logging.getLogger()
    dmin = shortest_distance(coord, cell) * tolerance
    pairs = pl.pairs_iter(coord, rc=dmin, cell=cell,
                          distance=False)
    for vtet, dtet in tetrahedra(pairs, dmin, coord, cell):
        p = dtet[0] + (dtet[1] + dtet[2] + dtet[3]) / 4
        p -= np.floor(p)
        yield p
