#!/usr/bin/env python3

#FAU Decoration of a 4-network
#読みこんだAR3Rの座標を、FAU構造における多面体vertexの位置とみなし、
#それらを連結するネットワークを六角柱で修飾して大きなネットワークを作る。

from math import acos, pi, sin, cos
from collections import defaultdict
import numpy as np
import logging
import re


def tune_angles(sixvecs, pivot):
    """
    Find the best origin of angles to make sum cos(3 th) largest
    """
    sixangles = []
    for i in range(len(sixvecs)):
        vec = sixvecs[i]
        cosine = np.dot(sixvecs[0],vec)
        if cosine > 1.0:
            cosine = 1.0
        angle = acos(cosine)
        sine   = np.cross(sixvecs[0], vec)
        if np.dot(sine, pivot) < 0:
            angle = -angle
        sixangles.append(angle)
    offset = 0
    while True:
        sum = 0.0
        dsum = 0.0
        for a in sixangles:
            sum += cos((a+offset)*6)
            dsum += -sin((a+offset)*6)
        doffset = dsum / 20.0
        if abs(doffset) < 1e-6:
            return offset
        offset += doffset
    

    
class decorate():
    def __init__(self, atoms, cell, pairs, Ncyl):
        """
        Ncyl is the number of cylinders to be inserted (>0)
        """
        #make netghbor list
        nei = defaultdict(set)
        for i,j in pairs:
            nei[i].add(j)
            nei[j].add(i)
        self.nei = nei
        self.atoms = atoms
        self.cell  = cell
        self.Ncyl  = Ncyl
        self.vertices = []
        self.fixedEdges = []
        for pair in pairs:
            self.one(pair)

    def one(self, pair):
        i,j = pair
        dij = self.atoms[j] - self.atoms[i]
        dij -= np.floor(dij + 0.5)
        dij = np.dot(dij, self.cell)
        scale = np.linalg.norm(dij) 
        dij /= scale
        sixpairs = []
        for k in self.nei[i]:
            if k != j:
                sixpairs.append((k,i))
        for k in self.nei[j]:
            if k != i:
                sixpairs.append((k,j))
        #Regularize the dihedral angles
        #to point them 6-fold directions.
        #by adding an offset
        sixvecs = []
        for j,k in sixpairs:
            vec = self.atoms[j] - self.atoms[k]
            vec -= np.floor(vec + 0.5)
            vec = np.dot(vec, self.cell)
            #orthogonalize
            shadow = np.dot(dij, vec)
            vec -= shadow*dij
            vec /= np.linalg.norm(vec)
            #print(np.dot(vec,dij))
            sixvecs.append(vec)
        offset = tune_angles(sixvecs, dij)
        offset += pi/6 #30 degree
        x = sixvecs[0].copy()
        z = dij
        y = np.cross(z,x)
        for j in range(6):
            a = j*pi*2/6 + offset
            sixvecs[j] = x*cos(a) + y*sin(a)
        #determine r
        #assume edge length is 1
        #the radius of the outer sphere of the polyhed is sqrt(3/2)
        L = (3/2)**0.5 * 2 + self.Ncyl
        r = 1/L     #edge len = radius of cyl
        rp = (3/2)**0.5 / L  # = radius of polyhed
        #
        a = np.dot(self.atoms[i],self.cell)
        s = ""
        for j in range(0, self.Ncyl+1):
            vec0 = dij*(rp + j*r)*scale + a
            for vec in sixvecs:
                rpos = vec0 + vec*r*scale
                pos = np.dot(rpos, np.linalg.inv(self.cell))
                self.vertices.append(pos)
            first = len(self.vertices)-6
            if j % 2 == 0:
                for k in range(5):
                    self.fixedEdges.append((first+k, first+k+1))
                self.fixedEdges.append((first+5,first))
            else:
                for k in range(5):
                    self.fixedEdges.append((first+k+1, first+k))
                self.fixedEdges.append((first,first+5))
            if j > 0:
                for k in range(6):
                    if k % 2 == 0:
                        self.fixedEdges.append((first+k, first+k-6))
                    else:
                        self.fixedEdges.append((first+k-6, first+k))


logger = logging.getLogger()
from genice.lattices import ice1c # base topology
cell1c = np.diag(np.fromstring(ice1c.cell, sep=" "))
waters1c = np.fromstring(ice1c.waters, sep=" ")
waters1c = waters1c.reshape((waters1c.shape[0]//3,3))
pairs1c = np.fromstring(ice1c.pairs, sep=" ", dtype=int)
pairs1c = pairs1c.reshape((pairs1c.shape[0]//2,2))


def argparser(arg):
    global Ncyl, coord, celltype, cell, waters, fixed
    assert re.match("^[0-9]+$", arg) is not None, "Argument must be an integer."
    Ncyl = int(arg)
    logger.info("Superlattice {0}xFAU".format(Ncyl))
    dec = decorate(waters1c, cell1c, pairs1c, Ncyl)
    coord='relative'
    celltype='rect'
    cell = "{0} {1} {2}".format(dec.cell[0,0],dec.cell[1,1],dec.cell[2,2])
    waters = dec.vertices
    fixed = dec.fixedEdges

# default.
argparser("1")
