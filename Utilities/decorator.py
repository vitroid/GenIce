#!/usr/bin/env python3

#FAU Decoration of a 4-network
#読みこんだAR3Rの座標を、FAU構造における多面体vertexの位置とみなし、
#それらを連結するネットワークを六角柱で修飾して大きなネットワークを作る。

from math import acos, pi, sin, cos
from collections import defaultdict
import sys
import itertools as it

import numpy as np

from genice import pairlist as pl
from genice import yaplotlib as yp

def LoadAR3R(file):
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if "@BOX3" == line[:5]:
            line  = file.readline()
            cols  = line.split()
            dimen = np.array([float(x) for x in cols])
            cell  = np.diag(dimen)
        elif "@AR3R" == line[:5]:
            natom = int(file.readline())
            Os = []
            for i in range(natom):
                line = file.readline()
                cols = line.split()
                xyz = np.array([float(x) for x in cols])
                xyz -= np.floor( xyz + 0.5 )
                Os.append(xyz)
            Os = np.array(Os)
    #the last line is cell shape
    return cell, Os

                 

def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99
    if pairs is None:
        iter = it.combinations(coord,2)
    else:
        iter = [(coord[i],coord[j]) for i,j in pairs]
    for c1,c2 in iter:
        d = c1-c2
        d -= np.floor(d + 0.5)
        r = np.dot(d,cell)
        rr = np.dot(r,r)
        if rr < dmin:
            dmin = rr
    return dmin**0.5


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
    vertices = []
    fixedEdges = []
    yap = ""
    def __init__(self, atoms, cell, pairs, Ncyl, mode = "python"):
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
        self.mode = mode
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
        s += yp.Color(3)
        s += yp.Line(a, a+dij*scale)
        for j in range(0, self.Ncyl+1):
            vec0 = dij*(rp + j*r)*scale + a
            for vec in sixvecs:
                s += yp.Line(vec0, vec0 + vec*r*scale)
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
        self.yap += s

    def __str__(self):
        if self.mode == "yaplot":
            return self._str_yaplot()
        else:
            return self._str_python()

    def _str_python(self):
        s = []
        s.append("coord='relative'")
        s.append("celltype='rect'")
        s.append("cell='{0} {1} {2}'".format(self.cell[0,0],self.cell[1,1],self.cell[2,2]))
        s.append('waters="""')
        for a in self.vertices:
            s.append("{0} {1} {2}".format(*a))
        s.append('"""')
        s.append('fixed="""')
        for i,j in self.fixedEdges:
            s.append("{0} {1}".format(i,j))
        s.append('"""')
        return "\n".join(s)


    def _str_yaplot(self):
        return self.yap

    

def main():
    mode = "python"  #output python module for GenIce
    if sys.argv[1] == "-y":
        mode = "yaplot"
        sys.argv.pop(1)
    Ncyl = int(sys.argv[1])
    #read O positions
    cell, atoms = LoadAR3R(sys.stdin)
    #build the network
    dmin = shortest_distance(atoms, cell)
    pairs = [v for v in pl.pairlist_crude(atoms, dmin*1.1, cell, distance=False)]

    #decorate
    dec = decorate(atoms, cell, pairs, Ncyl, mode)
    print(dec)

if __name__ == "__main__":
    main()
