from genice.libgenice import GenIce
from genice import openscad2
import numpy as np

class Formatter(GenIce):
    """
    cell is in nm

    openscad2 comes up with OO style
    """
    def run(self, options):
        water_type    = options.water[0]
        guests        = options.guests
        for d in range(3):
            self.rep[d] += 2  #Extend the size,then cut off later.

        self.stage1()   #replicate the unit cell
        self.stage2()   #prepare random graph

        scale=50
        roxy=0.07
        rbond=0.06
        self.logger.info("Output water molecules in OpenSCAD format revised.")
        rep = np.array(self.rep)
        trimbox    = self.cell *np.array([(rep[i]-2)/rep[i] for i in range(3)])
        trimoffset = (self.cell[0]+self.cell[1]+self.cell[2])/rep

        margin = 0.2 # expansion relative to the cell size
        lower = (1.0 - margin) / rep
        upper = (rep - 1.0 + margin) / rep
    
        bonds = []
        for i,j in self.graph.edges_iter(data=False):
            s1 =self.reppositions[i]
            s2 =self.reppositions[j]
            d = s2-s1
            d -= np.floor( d + 0.5 )
            self.logger.debug("Len {0}-{1}={2}".format(i,j,np.linalg.norm(d)))
            s2 = s1 + d
            if ( (lower[0] < s1[0] < upper[0] and lower[1] < s1[1] < upper[1] and lower[2] < s1[2] < upper[2] ) or
                (lower[0] < s2[0] < upper[0] and lower[1] < s2[1] < upper[1] and lower[2] < s2[2] < upper[2] ) ):
                bonds.append( (np.dot(s1,self.cell), np.dot(s2,self.cell)))
        
        nodes = []
        for s1 in self.reppositions:
            if lower[0] < s1[0] < upper[0] and lower[1] < s1[1] < upper[1] and lower[2] < s1[2] < upper[2]:
                nodes.append( np.dot(s1, self.cell) )
        
        o = openscad2.OpenScad()
        objs = [o.sphere(r="Roxy").translate(node) for node in nodes] + [o.bond(s1,s2,r="Rbond") for s1,s2 in bonds]
        #operations
        ops = [openscad2.bondfunc,
            o.defvar("$fn", 20),
            o.defvar("Roxy", roxy),
            o.defvar("Rbond", rbond),
            ( o.rhomb(trimbox).translate(trimoffset) & o.add(*objs) ).translate(-trimoffset).scale([scale,scale,scale])]
        s = o.encode(*ops)
        print(s,end="")
