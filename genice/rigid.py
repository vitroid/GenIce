#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
from math import *
import sys
import numpy as np
import numpy.linalg
import logging

#Rigid.py uses vertical vectors while other pars of genice uses horizontal vectors.

#v1,v2,v3: numpy arrays
#return value: fourth vector candidate
def decide1(v1,v2,v3):
    v4 = -(v1+v2+v3)
    return v4 / np.linalg.norm(v4)


def decide2(v1,v2):
    d = v1+v2
    d /= np.linalg.norm(d)
    v = np.cross(v1,v2)
    v /= np.linalg.norm(v)
    d *= -sqrt(1./3)
    v *=  sqrt(2./3)
    return d+v,d-v


#outer product of two vectors; return None if vector is too small
def op(i,j,check=True):
    if check and (np.linalg.norm(i) < 0.1 or np.linalg.norm(j) < 0.1):
        return np.array([])
    a = np.cross(i,j)
    if check and np.linalg.norm(a) < 0.1:
        return np.array([])
    return a



#calculate quaternions from a rotation matrix (three orthogonal vectors)
def rotmat2quat0(i,j,k):
    #print sqlen(i),sqlen(j),sqlen(k)
    ex = np.array((1.0, 0.0, 0.0))
    ey = np.array((0.0, 1.0, 0.0))
    ez = np.array((0.0, 0.0, 1.0))
    
    # i軸をx軸に移す回転の軸は、iとxの2分面上にある。
    # j軸をy軸に移す回転の軸は、jとyの2分面上にある。
    # そして、それらを同時にみたす回転の軸は、2つの2分面の交線である。
    # 交線は、2つの面の法線のいずれとも直交する=外積である。*/
    
    a = op(i-ex, j-ey)
    if a.size == 0:
        a = op(i-ex, k-ez)
        if a.size == 0:
            a = op(k-ez, j-ey)
            if a.size == 0:
                #sys.stderr.write("outer prod warning\n")
                #//全く回転しないケース
                return 1.0, 0.0, 0.0, 0.0
    a /= np.linalg.norm(a)
    #/*回転軸aが求まったので、x軸をi軸に重ねる回転の大きさを求める。。*/
    x0 = ex - a[0]*a
    i0 = i  - a[0]*a
    if np.linalg.norm(i0) < 0.1:
        i0 = j - a[1]*a
        x0 = ey - a[1]*a

    i0 /= np.linalg.norm(i0)
    x0 /= np.linalg.norm(x0)
    cosine = np.dot(i0,x0.transpose())
    if cosine < -1.0 or cosine > 1.0:
        cosh = 0.0
        sinh = 1.0
    else:
        cosh=sqrt((1.0+cosine)*0.5)
        sinh=sqrt(1.0-cosh*cosh)
    #//fprintf(stderr,"sinh %24.17e %24.17e %24.17e\n",cosine,cosh,sinh);
    #/*outer product to determine direction*/
    o = op(i0,x0,False)
    a = a.transpose()
    if np.dot(o,a) < 0.0:
        sinh=-sinh;
    
    return np.array((cosh, -sinh*a[0], +sinh*a[1], -sinh*a[2]))



def rotmat2quat(m):
    #print "rotmat2quat is not reliable yet."
    #sys.exit(1)
    n = m.transpose()
    return rotmat2quat0(np.array(n[0]),np.array(n[1]),np.array(n[2]))


def quat2rotmat(q):
    a,b,c,d = q
    sp11=(a*a+b*b-(c*c+d*d))
    sp12=-2.0*(a*d+b*c)
    sp13=2.0*(b*d-a*c)
    sp21=2.0*(a*d-b*c)
    sp22=a*a+c*c-(b*b+d*d)
    sp23=-2.0*(a*b+c*d)
    sp31=2.0*(a*c+b*d)
    sp32=2.0*(a*b-c*d)
    sp33=a*a+d*d-(b*b+c*c)
    return np.array([[sp11,sp12,sp13], [sp21,sp22,sp23], [sp31,sp32,sp33]])



def euler2quat(e):
    ea,eb,ec = e
    a=cos(ea/2)*cos((ec+eb)/2)
    b=sin(ea/2)*cos((ec-eb)/2)
    c=sin(ea/2)*sin((ec-eb)/2)
    d=cos(ea/2)*sin((ec+eb)/2)
    return np.array((a,b,c,d))
    


def quat2euler(q):
    e = np.zeros(3)
    if q[0] == 1.0:
        return e
    p = 2*(q[0]**2+q[3]**2) - 1
    if p>1.0:
        p=1.0
    if p<-1.0:
        p=-1.0
    e[0] = acos(p)
    thh = e[0] / 2.0
    sinthh = sin(thh)
    costhh = cos(thh)
    p = q[0]/costhh
    if p>1.0:
        p=1.0
    if p<-1.0:
        p=-1.0
    p = acos(p)
    if sinthh == 0.0:
        s = 1.0
    else:
        s = q[1]/sinthh
    if s>1.0:
        s=1.0
    if s<-1.0:
        s=-1.0
    s=acos(s)
    if q[3]<0.0:
        p=2*pi-p
    if q[2]>0:
        e[2] = p+s
        e[1] = p-s
    else:
        e[2] = p-s
        e[1] = p+s
    e[1:3] %= (2*pi)
    return e


def qadd(q1,q2):
    a1,b1,c1,d1 = q1
    a2,b2,c2,d2 = q2
    a3=a1*a2-b1*b2-c1*c2-d1*d2
    b3=a1*b2+b1*a2+c1*d2-d1*c2
    c3=a1*c2+c1*a2-b1*d2+d1*b2
    d3=a1*d2+d1*a2+b1*c2-c1*b2
    if a3 < 0:
       a3 = -a3
       b3 = -b3
       c3 = -c3
       d3 = -d3
    return np.array([a3,b3,c3,d3])



def qmul(q,x):
    if q[0] >= 1.0:
        return np.array(q)
    phi = acos(q[0])
    sine = sqrt(1.0-q[0]**2)
    if q[1] < 0:
        phi=-phi
        sine=-sine
    phi *= x
    sine = sin(phi) / sine
    a = np.zeros(4)
    a[0] = cos(phi)
    a[1] = q[1]*sine
    a[2] = q[2]*sine
    a[3] = q[3]*sine
    return a


def test_rotation():
    e = np.array((0.2,0.3,0.4))
    print(e)
    print(quat2euler(euler2quat(e)))
    print(euler2quat(e))
    print(rotmat2quat(quat2rotmat(euler2quat(e))))
    q = euler2quat(e)
    q = qmul(q, 0.5)
    print(qadd(q,q))


    

def test():
    """Testing Docstring"""
    q1 = np.array([0.5,0.5,0.5,0.5]);
    q2 = np.array([0.5,0.5,0.5,-0.5]);
    q21 = qadd(q2,qmul(q1,-1))
    print(q21)
    print(qadd(q21,q1))
    q3 = np.array([0.60876143, -0.45804276,  0.45804276, -0.45804276])
    print(qadd(q3,q1))
    print(qadd(q1,q3))
    test_rotation()
    #first molecule by genice 1h -r 1 1 1 --format q
    q4 = np.array([0.0692,-0.2480,0.2911,-0.9214])
    e4 = quat2euler(q4)
    print(e4)
    qe4 = euler2quat(e4)
    print(qe4)


#http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
def rand_rotation_matrix(deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = np.random.uniform(size=(3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
        np.sin(phi) * r,
        np.cos(phi) * r,
        np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M



if __name__ == "__main__":
    test()
