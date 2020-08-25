# coding: utf-8

desc={"ref": {"Codes": "https://github.com/vitroid/Yaplot"},
      "brief": "Yaplot.",
      "usage": """
Usage: genice icename -f yaplot[options]

options:
    H=x   Set the radius of H to be x.
"""}



from collections import defaultdict
from logging import getLogger
import numpy as np
import yaplotlib as yp


import genice.formats
class Format(genice.formats.Format):
    size_H = 0.01


    def __init__(self, **kwargs):
        unknown = dict()
        for k, v in kwargs.items():
            if k == "H":
                self.size_H = float(v)
            else:
                unknown[k] = v
        super().__init__(**unknown)


    def hooks(self):
        return {1:self.hook1,
                2:self.hook2,
                6:self.hook6,
                7:self.hook7}


    def hook1(self, ice):
        logger = getLogger()
        logger.info("Hook1: Draw the cell in Yaplot format.")
        s = yp.Layer(2)
        x,y,z = ice.repcell.mat
        for p,q,r in ((x,y,z),(y,z,x),(z,x,y)):
            for a in (np.zeros(3), p, q, p+q):
                s += yp.Line(a,a+r)
        print(s,end="")
        logger.info("Hook1: end.")


    def hook2(self, ice):
        logger = getLogger()
        if self.size_H > 0:
            return
        logger.info("Hook2: Output CoM of water molecules in Yaplot format.")
        # prepare the reverse dict
        waters = defaultdict(dict)
        pos = ice.reppositions @ ice.repcell.mat
        s = ""
        for p in pos:
            s += yp.Layer(4)
            s += yp.Color(3)
            s += yp.Size(0.03)
            s += yp.Circle(p)
        print(s, end="")
        logger.info("Hook2: end.")
        return True


    def hook6(self, ice):
        logger = getLogger()
        logger.info("Hook6: Output water molecules in Yaplot format.")
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        # prepare the reverse dict
        waters = defaultdict(dict)
        for atom in ice.atoms:
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
            s += yp.Size(self.size_H)
            s += yp.Circle(H0)
            s += yp.Circle(H1)
            s += yp.Color(2)
            s += yp.Layer(1)
            s += yp.Text(O, "{0}".format(order))
        s += yp.Layer(3)
        s += yp.Color(4)
        s += yp.ArrowType(1)
        s += yp.Size(0.03)
        for i,j in ice.spacegraph.edges(data=False):
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
        self.nwateratoms = len(ice.atoms)
        logger.info("Hook6: end.")


    def hook7(self, ice):
        logger = getLogger()
        logger.info("Hook7: Output water molecules in Yaplot format.")
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        gatoms = ice.atoms[self.nwateratoms:]
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
        s = '#' + "\n#".join(ice.doc) + "\n" + s
        print(s)
        logger.info("Hook7: end.")
