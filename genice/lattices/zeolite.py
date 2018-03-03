#!/usr/bin/env python3

#system modules
import os
import sys
import itertools as it
import logging
#external modules
import numpy as np
from requests import get #requests package
import validators        #validators package
from cif2ice import read_cif

def shortest_distance(atoms):
    dmin = 1e99
    for a1,a2 in it.combinations(atoms,2):
        name1,x1,y1,z1 = a1
        name2,x2,y2,z2 = a2
        d = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        if d < dmin:
            dmin = d
    return dmin**0.5


def is_unique(L, pos):
    for x in L:
        d = x-pos
        d -= np.floor(d + 0.5)
        if np.dot(d,d) < 0.0000001:
            return False
    return True


def genice_lattice(atoms, box, matchfunc=None):
    global cell, celltype, waters, coord, density
    logger = logging.getLogger()
    filtered = []
    if matchfunc is not None:
        for a in atoms:
            if matchfunc(a[0]):
                filtered.append(a)
    else:
        filtered = atoms
    dmin = shortest_distance(filtered)
    scale = 2.76 / dmin

    celltype = "triclinic"
    if (len(box) == 6):
        cell = np.array([[box[0],0,0],[box[1],box[2],0],[box[3],box[4],box[5]]])
    else:
        cell = np.diag(box)
    volume = np.linalg.det(cell)
    icell  = np.linalg.inv(cell)
    uniques = []
    for name,x,y,z in filtered:
        rpos = np.dot([x,y,z], icell)
        rpos -= np.floor(rpos)
        #Do twice to reduce the floating point uncertainty.
        #(Hint: assume the case when x=-1e-33.)
        rpos -= np.floor(rpos)
        if is_unique(uniques, rpos):
            uniques.append(rpos)
    waters = uniques
    coord = "relative"
    # bondlen = 3
    density = len(filtered)*18.0/(volume*scale**3*1e-24*6.022e23)

    
def download(url, file_name):
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)
        
def argparser(arg):
    logger = logging.getLogger()
    args = arg.split(":")
    Nbox = (1,1,1) # should be given as an option
    name = args[0]
    logger.info(__name__)
    #input must be a file......too bad.
    if os.path.exists(name):
        fNameIn = name
    else:
        if validators.url(name):
            URL = name
            name = os.path.basename(name)
            if name[-4:] in (".cif", ".CIF"):
                name = name[:-4]
        elif __name__[-16:] == "lattices.zeolite":
            # it only works when my module name is zeolite.
            URL = "http://www.iza-structure.org/IZA-SC/cif/"+name+".cif"
        fNameIn = name + ".cif"
        assert validators.url(URL)
        download(URL, fNameIn)
    logger.info("Input: {0}".format(fNameIn))
    atoms, box = read_cif.read_and_process(fNameIn, Nbox, make_rect_box=False)
    genice_lattice(atoms, box, matchfunc=lambda x: x[0] != "O")


        
# Do nothing by default; it causes an error when arguments are missing.
