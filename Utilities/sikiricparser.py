#!/usr/bin/env python3

import sys
import numpy as np
from math import pi, sin, cos, sqrt,floor
import logging
import itertools as it
import triangles as ta
from collections import defaultdict

def is_zero(x):
    return -1e-14 < x < 1e-14


def blurred(v, wrapped=True):
    yield v  #in advance
    w=2
    for x in range(-w,w+1):
        for y in range(-w,w+1):
            for z in range(-w,w+1):
                if wrapped:
                    yield (v[0]+x)%100000, (v[1]+y)%100000, (v[2]+z)%100000
                else:
                    yield v[0]+x, v[1]+y, v[2]+z


#L must be a wrapped list.
def find_with_error(v, L, wrapped=True):
    #logger = logging.getLogger()
    #debug = (v == (-17176.0, 17176.0, 6001.0))
    for b in blurred(v, wrapped):
        #if debug:
        #    logger.debug((b, b in L))
        if b in L:
            return b
    return None


encoded = []

#for water. water is base 400000
def encode(r):  #base 400000
    logger= logging.getLogger()
#    r /= 800.0
    
#    x = (floor(r[0]+0.5) + 500) % 500 *800
#    y = (floor(r[1]+0.5) + 500) % 500 *800
#    z = (floor(r[2]+0.5) + 500) % 500 *800
    r /= 4
    x = (floor(r[0]+0.5) + 100000) % 100000. * 4
    y = (floor(r[1]+0.5) + 100000) % 100000. * 4
    z = (floor(r[2]+0.5) + 100000) % 100000. * 4
    p = np.array([x,y,z])
    dmin = 1e99
    qmin = 0
    for q in encoded:
        d = p - q
        d -= np.floor( d / 400000 + 0.5 ) * 400000
        Ld = np.dot(d,d)
        if Ld < dmin:
            dmin = Ld
            qmin = q
    if dmin > 100:
        encoded.append(p)
        return (x,y,z)
    return tuple(qmin)

#for water
def shortest_distance(relpos, box):
    logger = logging.getLogger()
    dmin = 1e99
    amin = []
    for a1,a2 in it.combinations(relpos,2):
        d = a1-a2
        d -= np.floor( d + 0.5 )
        d = np.dot(d, box)
        dd = np.dot(d,d)
        if dd < dmin:
            dmin = dd
            amin = a1,a2
    logger.debug("amin:{0}".format(amin))
    return dmin**0.5



def lineparser(line):
    columns = line.split()
    el = columns[0].split("=")
    mode = el[0]
    if mode == "cellstructure":
        return mode,0,[float(x) for x in [el[1],]+columns[1:6]]
    order = int(el[1])
    #print line
    if columns[2] == "0,":
        xnum = 0
    elif columns[2] == "1,":
        xnum = 100000
    elif columns[2] == "-1,":
        xnum = -100000
    elif columns[2] == "2,":
        xnum = 200000
    else:
        xnum, xden = columns[2].split("/")
        xden = xden[:-1]
        xnum = int(xnum) * 100000 // int(xden)
    if columns[3] == "0,":
        ynum = 0
    elif columns[3] == "1,":
        ynum = 100000
    elif columns[3] == "-1,":
        ynum = -100000
    elif columns[3] == "2,":
        ynum = 200000
    else:
        ynum, yden = columns[3].split("/")
        yden = yden[0:-1]
        ynum = int(ynum) * 100000 // int(yden)
    if columns[4] == "0":
        znum = 0
    elif columns[4] == "1":
        znum = 100000
    elif columns[4] == "-1":
        znum = -100000
    elif columns[4] == "2":
        znum = 200000
    else:
        znum, zden = columns[4].split("/")
        znum = int(znum) * 100000 // int(zden)
    #xnum %= 100000
    #ynum %= 100000
    #znum %= 100000
    return mode,order-1,np.array([float(v) for v in (xnum,ynum,znum)])

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s %(levelname)s %(message)s")
    
    logger = logging.getLogger()
    x0=0
    y0=0
    z0=0
    o0=0
    bonddest = dict()
    box = []
    atoms = dict()
    for line in sys.stdin:
        mode,order,r = lineparser(line)
        if mode == "iVect":
            r0= r
            rp = r%100000
            atoms[tuple(rp)] = order
            o0 = order
            bonddest[order] = dict()
        elif mode == "cellstructure":
            box = r
            #print r
        else:
            #iAdj
            r -= r0
            bonddest[o0][tuple(r)] = -1

    #Make the neighbor dict.
    #for apos in atoms:
    #    aord = atoms[apos]
    #    logger.debug("{0} pos {1}".format(aord,apos))
    for apos in atoms:
        aord = atoms[apos]
        for bvec in bonddest[aord]:
            cpos = np.array(apos) + np.array(bvec)
            cpos = cpos % 100000
            correct = find_with_error(tuple(cpos), atoms)
            if correct is not None:
                cord = atoms[tuple(correct)]
                bonddest[aord][bvec] = cord
            else:
                #test if there is really no vertex.
                dmin = 1e99
                dv   = None
                c2 = cpos
                for a2 in atoms:
                    a3 = np.array(a2)
                    d  = a3-c2
                    dd = np.dot(d,d)
                    if dd < dmin:
                        dmin = dd
                        dv = d
                logger.debug("Closest:{0} {1} Center[{5}] {2} + delta {3} = cpos {4}".format(dmin,dv,apos,bvec,cpos,aord))
                
    
    alpha = box[3]*pi/180
    beta  = box[4]*pi/180
    gamma = box[5]*pi/180
    A = np.array([1.0, 0.0, 0.0])
    Rg = np.array([[cos(gamma),  sin(gamma), 0],
                   [-sin(gamma), cos(gamma), 0],
                   [0,0,1]])
    B = np.dot(A,Rg)
    #print(B,box[5])

    #unknown C
    #A.C = cos(beta)
    #B.C = cos(alpha)
    #C.C = 1
    Cx = cos(beta)
    #Bx*Cx+By*Cy=cos(alpha)
    Cy = (cos(alpha)-B[0]*Cx)/B[1]
    Cz = sqrt(1-Cx**2-Cy**2)
    C = np.array([Cx,Cy,Cz])
    #print(A,B,C)
    #print(np.dot(B,C),cos(alpha))
    #print(np.dot(C,A),cos(beta ))
    #print(np.dot(A,B),cos(gamma))
    #sys.exit(1)
    A *= box[0]
    B *= box[1]
    C *= box[2]
    cell = np.array([A,B,C])

    s = ""
    #for key in atoms:
    #    logger.debug("[{0},{1},{2}]".format(key[0], key[1], key[2]))
        
    cagecounts = dict()
    for i in (12,14,15,16):
        cagecounts[i] = 0
    for a in atoms:
        na= np.array(a)
        i = atoms[a]
        nnei = len(bonddest[i])
        cagecounts[nnei] += 1
    logger.debug("cagecounts:{0}".format(cagecounts))
    #for a in atoms:
    #    na   = np.array(a)   #location of the atom
    #    aord = atoms[a]      #label of the atom
    #    #Test if vectors are isotropically distributed
    #    sum = np.zeros(3)
    #    for bvec in bonddest[aord]:
    #        sum += np.array(bvec)
    #    logger.debug("Balance {0}".format(sum / len(bonddest[aord])))
    dones = set()
    #done2 = set()
    dualedges = set()
    waters = set()
    nv = 0
    for a in atoms:
        na   = np.array(a)   #location of the atom
        aord = atoms[a]      #label of the atom
        assert len(bonddest[aord]) in (12,14,15,16)   #FK rule
        #Prepare the adjacent pair list
        nei = defaultdict(list)
        #check if there is the neighbor vertex
        for bvec in bonddest[aord]:
            cpos = na + np.array(bvec)
            cpos %= 100000
            assert find_with_error(tuple(cpos), atoms), bvec
        #Find the neighboring vector pairs originating from aord
        nedge = 0
        for bvec in bonddest[aord]:
            bord = bonddest[aord][bvec]
            #logger.debug("B nei:{0}".format(len(bonddest[bord])))
            nedge1 = 0
            #logger.debug(bonddest[aord])
            for cvec in bonddest[bord]:
                cpos = tuple(np.array(bvec) + np.array(cvec))
                correct = find_with_error(cpos, bonddest[aord], wrapped=False)
                #logger.debug(("<>",cpos, correct))
                if correct is not None:
                    #both correct and bvec are in bonddest[aord]
                    #There is a triangle
                    nei[bvec].append(correct)
                    nei[correct].append(bvec)
                    nedge +=1
                    nedge1 += 1
                                
                    #test if there is really no vertex.
                    #dmin = 1e99
                    #c2 = np.array(cpos)
                    #for b2 in bonddest[aord]:
                    #    b3 = np.array(b2)
                    #    d  = b3-c2
                    #    dd = np.dot(d,d)
                    #    if dd < dmin:
                    #        dmin = dd
                    #logger.debug("Closest:{0}".format(dd))
                    #pass
            #logger.debug("# of edges at {1}: {0}".format(nedge1, bvec))
        #logger.debug("NVertex:{1} NEdge:{0}".format(nedge, len(bonddest[aord])))
        #for bvec in nei:
        #    logger.debug(nei[bvec])
        #triangles sharing the edge
        edgesharers = defaultdict(list)
        for tri in ta.triangles(nei):
            j,k,l = tri
            #logger.debug(tri)
            edgesharers[frozenset([j,k])].append(set(tri))
            edgesharers[frozenset([k,l])].append(set(tri))
            edgesharers[frozenset([l,j])].append(set(tri))
        waters0 = set()
        dualedges0 = set()
        for edge in edgesharers:
            assert len(edgesharers[edge]) == 2, "The edge is shared by {0} trianglular faces.".format(len(edgesharers[edge]))
            j,k = edge
            l,m = (edgesharers[edge][0] | edgesharers[edge][1]) - edge
            #triangles jkl and jkm are adjacent.
            c1 = na*4 + np.array(j) + np.array(k) + np.array(l)
            c2 = na*4 + np.array(j) + np.array(k) + np.array(m)
            ic1 = encode(c1)
            ic2 = encode(c2)
            dualedges0.add(frozenset((ic1,ic2)))
            waters0.add(ic1)
            waters0.add(ic2)
        #logger.debug("vertices: {0} waters: {1} edges: {2}".format(len(bonddest[aord]),len(waters0),len(dualedges0)))
        #Euler characteristics
        assert len(bonddest[aord])*2 - len(waters0) == 4
        assert (len(bonddest[aord]) - 12)*6 + 12*5 == len(dualedges0)*2
        waters |= waters0
        dualedges |= dualedges0


    waterids = dict()
    count = 0
    for water in waters:
        waterids[water] = count
        count += 1
    s += 'pairs="""\n'
    for edge in dualedges:
        i,j = edge
        s += '{0} {1}\n'.format(waterids[i],waterids[j])
    s += '"""\n\n'
    #nv = len(atoms)
    #assert len(waters)*4 == nv, "{0} {1}".format(len(waters)*4, nv)
    for water in sorted(waters):
        logger.debug(water)
    waters = np.array([water for water in waters]) / 400000
    #for water in sorted(waters):
    #    print(water)
    dmin = shortest_distance(waters, cell)
    logger.debug("dmin:{0}".format(dmin))
    assert dmin > 0.1
    scale = 2.76 / dmin
    s += 'waters="""\n'
    for water in waters:
        s += '{0} {1} {2}\n'.format(*water)
    s += '"""\n\n'
    s += 'coord= "relative"\n\n'
    s += 'cages="""\n'
    for a in atoms:
        na= np.array(a)
        i = atoms[a]
        s += '{0} {1} {2} {3}\n'.format(len(bonddest[i]), *(np.array(a)/100000))
    s += '"""\n\n'
    s += "bondlen = 3\n\n"
    cell *= scale
    A,B,C = cell
    if is_zero(A[1]) and is_zero(A[2]) and is_zero(B[0]) and is_zero(B[2]) and is_zero(C[0]) and is_zero(C[1]):
        s += "celltype = 'rect'\n\n"
        s += 'cell = """\n'
        s += '{0} {1} {2}\n'.format(A[0],B[1],C[2])
        s += '"""\n\n'
    else:
        s += "celltype = 'triclinic'\n\n"
        s += 'cell = """\n'
        s += '{0} {1} {2}\n'.format(*A)
        s += '{0} {1} {2}\n'.format(*B)
        s += '{0} {1} {2}\n'.format(*C)
        s += '"""\n\n'
    volume = np.linalg.det(cell)
    density = len(waters)*18.0/(volume*1e-24*6.022e23)
    s += "density = {0}\n\n".format(density)
    header = '"""\n'
    header += "Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.\n\n"
    header += "Cage composition:\n (12,14,15,16) = ("
    for c in (12,14,15,16):
        header += "{0},".format(cagecounts[c])
    header += ")\n"
    header += '"""\n\n'
    print(header+s)
    
if __name__ == "__main__":
    main()
