# coding: utf-8
"""
Povray format.
"""

from collections import defaultdict
import numpy as np

from genice import rigid


def Block(name, content):
    return " " + name + " { " + content + " } "

def Vector(v):
    return "<{0:.3f}, {1:.3f}, {2:.3f}> ".format(*v)

def Juxtapose(v):
    return ",".join(v)

def Atom(atomtype, pos):
    return Block( "sphere", Juxtapose( [Vector(pos), "R{0}".format(atomtype) ] ) + Block( "material", "MAT{0}".format(atomtype) ) ) + "\n"
    
def Bond(bondtype, pos1, pos2):
    return Block( "cylinder", Juxtapose( [Vector(pos1), Vector(pos2), "R{0}".format(bondtype)] ) + Block( "material", "MAT{0}".format(bondtype) ) ) + "\n"
    

def hook6(lattice):
    global nwateratoms
    lattice.logger.info("Hook6: Output water molecules in Povray format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    # prepare the reverse dict
    waters = defaultdict(dict)
    for atom in lattice.atoms:
        resno, resname, atomname, position, order = atom
        if "O" in atomname:
            waters[order]["O"] = position
        elif "H" in atomname:
            if "H0" not in waters[order]:
                waters[order]["H0"] = position
            else:
                waters[order]["H1"] = position
    s = '#include "default.inc"\n'
    s += "union {\n"
    for order, water in waters.items():
        O = water["O"]
        H0 = water["H0"]        
        H1 = water["H1"]        
        s += Atom("O",O)
        s += Atom("H",H0)
        s += Atom("H",H1)
        s += Bond("OH",O,H0)
        s += Bond("OH",O,H1)
    for i,j in lattice.graph.edges_iter(data=False):
        if i in waters and j in waters:  # edge may connect to the dopant
            O = waters[j]["O"]
            H0 = waters[i]["H0"]
            H1 = waters[i]["H1"]
            d0 = H0 - O
            d1 = H1 - O
            rr0 = np.dot(d0,d0)
            rr1 = np.dot(d1,d1)
            if rr0 < rr1 and rr0 < 0.245**2:
                s += Bond("HB",H0,O)
            if rr1 < rr0 and rr1 < 0.245**2:
                s += Bond("HB",H1,O)
    print(s, end="")
    nwateratoms = len(lattice.atoms)
    lattice.logger.info("Hook6: end.")


def hook7(lattice):
    global nwateratoms
    lattice.logger.info("Hook7: Output water molecules in Povray format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    gatoms = lattice.atoms[nwateratoms:]
    s = ""
    H = []
    O  = ""
    for atom in gatoms:
        resno, resname, atomname, position, order = atom
        s += Atom(atomname,position)
    s = '//' + "\n//".join(lattice.doc) + "\n" + s
    s += "  translate " + Vector( -(lattice.repcell[0,:]+lattice.repcell[1,:]+lattice.repcell[2,:])/2 ) + "\n}\n"
    print(s)
    lattice.logger.info("Hook7: end.")
    

hooks = {7:hook7, 6:hook6}
