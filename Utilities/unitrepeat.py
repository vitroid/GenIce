#!/usr/bin.env python3

import sys
sys.path.append('../lib')
sys.path.append('../Prepare')
import fileloaders as fl
import numpy as np
#import yaplotlib as yp
import itertools as it
MAXCELL=11

def OccupyCell(nv, cells, ijk):
    if nv not in cells:
        #print("kill:",nv)
        cells.add(nv)
        for d in it.product(range(-1,2),range(-1,2),range(-1,2)):
            vv = tuple(np.dot(d,ijk)+nv)
            if 0 <= vv[0] < MAXCELL and 0 <= vv[1] < MAXCELL and 0 <= vv[2] < MAXCELL:
                #print("D:",d)
                OccupyCell(vv, cells, ijk)



def FindEmptyCell(nv, cells, box, r, rL, coord, ijk):
    ncell = 0
    if nv not in cells:
        ncell += 1
        for xyz in coord:
            xxv = xyz + nv*box
            p0v = np.dot(xxv,r)
            pv  = p0v / rL / rL
            pv -= np.rint(pv)
            prv = pv*rL
            print(prv[0],prv[1],prv[2])
        #print("NV:",nv)
        OccupyCell(nv, cells, ijk)
        xmin = max(nv[0]-1,0)
        xmax = min(nv[0]+2,MAXCELL)
        ymin = max(nv[1]-1,0)
        ymax = min(nv[1]+2,MAXCELL)
        zmin = max(nv[2]-1,0)
        zmax = min(nv[2]+2,MAXCELL)
        for xi,yi,zi in it.product(range(xmin,xmax),range(ymin,ymax),range(zmin,zmax)):
            ncell += FindEmptyCell((xi,yi,zi), cells, box, r, rL, coord, ijk)
    return ncell



iv = np.array([int(v) for v in sys.argv[1:4]])
jv = np.array([int(v) for v in sys.argv[4:7]])
kv = np.array([int(v) for v in sys.argv[7:10]])
yap = False
while True:
    line = sys.stdin.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if len(columns) > 0:
        if columns[0] =="@BOX3":
            box = fl.LoadBOX3(sys.stdin)
        elif columns[0] == "@AR3A":
            coord = fl.LoadAR3A(sys.stdin)

ii = np.dot(iv,iv)
jj = np.dot(jv,jv)
kk = np.dot(kv,kv)
ijk = np.array([iv,jv,kv])

riv = iv*box
ri  = np.linalg.norm(riv)
rih = ri * 0.5

rjv = jv*box
rj  = np.linalg.norm(rjv)
rjh = rj * 0.5

rkv = kv*box
rk  = np.linalg.norm(rkv)
rkh = rk * 0.5

r = np.array([ri,rj,rk])
rL = np.array([ri,rj,rk])

vv = np.cross(iv,jv)
vol = abs(np.dot(vv,kv))
print("Vol:",vol)
if yap:
    pass

cells = set()
ncell = FindEmptyCell((0,0,0), cells, box, r, rL, coord, ijk)
if vol != ncell:
    print("Internal error")
    sys.exit(1)
