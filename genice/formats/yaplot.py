# coding: utf-8

desc={"ref": {"Codes": "https://github.com/vitroid/Yaplot"},
      "brief": "Yaplot.",
      "usage": """
Usage: genice icename -f yaplot[options]

options:
    H=x   Set the radius of H to be x.
"""}



from collections import defaultdict
import numpy as np
import yaplotlib as yp


def hook0(lattice, args):
    lattice.logger.info("Hook0: ArgParser.")
    lattice.yaplot_size_H = 0.01
    if args != "":
        for arg in args.split(":"):
            cols = arg.split("=")
            if len(cols) == 2:
                if cols[0] == "H":
                    lattice.yaplot_size_H = float(cols[1])
                    continue
            assert False, "Unknown option."
    lattice.logger.info("Hook0: end.")
    
    



def hook1(lattice):
    lattice.logger.info("Hook1: Draw the cell in Yaplot format.")
    s = yp.Layer(2)
    x,y,z = lattice.repcell.mat
    for p,q,r in ((x,y,z),(y,z,x),(z,x,y)):
        for a in (np.zeros(3), p, q, p+q):
            s += yp.Line(a,a+r)
    print(s,end="")
    lattice.logger.info("Hook1: end.")


def hook2(lattice):
    global nwateratoms
    if lattice.yaplot_size_H > 0:
        return
    lattice.logger.info("Hook2: Output CoM of water molecules in Yaplot format.")
    # prepare the reverse dict
    waters = defaultdict(dict)
    pos = lattice.reppositions @ lattice.repcell.mat
    s = ""
    for p in pos:
        s += yp.Layer(4)
        s += yp.Color(3)
        s += yp.Size(0.03)
        s += yp.Circle(p)
    print(s, end="")
    lattice.logger.info("Hook2: end.")
    return True


def hook6(lattice):
    global nwateratoms
    lattice.logger.info("Hook6: Output water molecules in Yaplot format.")
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
    s = ""
    s += yp.Color(3)
    for order, water in waters.items():
        O = water["O"]        
        H0 = water["H0"]        
        H1 = water["H1"]        
        s += yp.Layer(4)
        s += yp.Color(3)
        s += yp.Size(0.03)
        s += yp.Circle(O)
        s += yp.Line(O,H0)
        s += yp.Line(O,H1)
        s += yp.Size(lattice.yaplot_size_H)
        s += yp.Circle(H0)
        s += yp.Circle(H1)
        s += yp.Color(2)
        s += yp.Layer(1)
        s += yp.Text(O, "{0}".format(order))
    s += yp.Layer(3)
    s += yp.Color(4)
    s += yp.ArrowType(1)
    s += yp.Size(0.03)
    for i,j in lattice.spacegraph.edges(data=False):
        if i in waters and j in waters:  # edge may connect to the dopant
            O = waters[j]["O"]
            H0 = waters[i]["H0"]
            H1 = waters[i]["H1"]
            d0 = H0 - O
            d1 = H1 - O
            rr0 = np.dot(d0,d0)
            rr1 = np.dot(d1,d1)
            if rr0 < rr1 and rr0 < 0.245**2:
                s += yp.Arrow(H0,O)
            if rr1 < rr0 and rr1 < 0.245**2:
                s += yp.Arrow(H1,O)
    print(s, end="")
    nwateratoms = len(lattice.atoms)
    lattice.logger.info("Hook6: end.")


def hook7(lattice):
    global nwateratoms
    lattice.logger.info("Hook7: Output water molecules in Yaplot format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    gatoms = lattice.atoms[nwateratoms:]
    palettes = dict()
    s = ""
    s += yp.Layer(4)
    s += yp.ArrowType(1)
    H = []
    O  = ""
    for atom in gatoms:
        resno, resname, atomname, position, order = atom
        if atomname in palettes:
            pal = palettes[atomname]
        else:
            pal = 4 + len(palettes)
            palettes[atomname] = pal
        s += yp.Color(pal)
        s += yp.Size(0.04)
        s += yp.Circle(position)
    s = '#' + "\n#".join(lattice.doc) + "\n" + s
    print(s)
    lattice.logger.info("Hook7: end.")
    

hooks = {0:hook0, 7:hook7, 6:hook6, 1:hook1, 2:hook2}
