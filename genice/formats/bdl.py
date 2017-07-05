# coding: utf-8
import numpy as np

def nearly_zero(x):
    return np.dot(x,x) < 1e-10


def run(lattice, water_type="TIP3P", guests=[]):
    """
    Scigress BDL file format (experimental)
    """
    lattice.stage1()   #replicate the unit cell
    lattice.stage2()   #prepare random graph
    lattice.stage3()   #Make an ice graph
    lattice.stage4()   #Depolarize
    lattice.stage5()   #Orientation
    lattice.stage6(water_type)  #Water atoms
    lattice.stage7B(guests)      #Guest atoms
    #lattice.logger.info("Total number of atoms: {0}".format(len(lattice.atoms)))
    lattice.logger.info("Output in Scigress BDL format.")

    a = lattice.cell[0,:]
    b = lattice.cell[1,:]
    c = lattice.cell[2,:]
    aL= np.linalg.norm(a)
    bL= np.linalg.norm(b)
    cL= np.linalg.norm(c)
    ab = np.dot(a,b)
    bc = np.dot(b,c)
    ca = np.dot(c,a)
    alpha = acos(bc/(bL*cL)) * 180 / pi
    beta  = acos(ca/(cL*aL)) * 180 / pi
    gamma = acos(ab/(aL*bL)) * 180 / pi
    automatic_mass = True
    s = []
    # Lines 1..2: Cell info
    s.append("{0:8.4f}  {1:8.4f}  {2:8.4f}".format(aL*10,bL*10,cL*10))
    s.append("{0:8.4f}  {1:8.4f}  {2:8.4f}".format(alpha,beta,gamma))
    #parse -g options again to determine the number of guest types
    guesttypes = set()
    for arg in guests:
        key, value = arg[0].split("=")
        guesttypes.add(value)
    # Line 3: number of molecular types
    s.append("{0:2}".format(len(guesttypes)+1))
    # Line 4: Guests.
    for gname in lattice.guestAtoms:
        guest = safe_import("molecule", gname)
        #molecule name, number of molecules, sites in a molecule, bonds in a molecule
        #last one is undefined i GenIce.
        s.append("{0:15}                {1:4}  {2:3}  {3:3}".format(gname, lattice.nGuestAtoms[gname], len(guest.sites), 0))
        for i in range(len(lattice.guestAtoms[gname])):
            molorder, resname, atomname, position = lattice.guestAtoms[gname][i]
            q = 0.0  #charge
            s.append("{0:4}  {1:8.4f}  {2:15.9f}  {3:15.9f}  {4:15.9f}".format(atomname,q,position[0],position[1],position[2]))

    #water
    water = safe_import("molecule", water_type)
    s.append("{0:15}                {1:4}  {2:3}  {3:3}".format(water_type, len(lattice.reppositions), len(water.sites), 0))
    for i in range(len(lattice.atoms)):
        molorder, resname, atomname, position = lattice.atoms[i]
        q = 0.0
        s.append("{0:4}  {1:8.4f}  {2:15.9f}  {3:15.9f}  {4:15.9f}".format(atomname,q,position[0],position[1],position[2]))

    for i in range(len(s)):
        line = s[i]
        print("{0:06}    ".format(i+1) + line)

#
#limitations:            
#1)It is a structure file and force field inf should not be included.
#however, BDL requests the partial carge info in it.
#2)Atomic naming rule is always a problem.
#
