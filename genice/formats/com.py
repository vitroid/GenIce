from genice.formats.baseclass import GenIce
from genice import rigid
import numpy     as np


class Formatter(GenIce):
    """
    Centers of mass of water molecule 
    """
    def run(self, options):
        water_type    = options.water[0]
        guests        = options.guests
        self.stage1()   #replicate the unit cell
        self.logger.info("Output centers of mass of water molecules.")
        s = ""
        if self.cell[1,0] == 0 and self.cell[2,0] == 0 and self.cell[2,1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(self.cell[0,0]*10,self.cell[1,1]*10,self.cell[2,2]*10)
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(self.cell[0,d]*10,self.cell[1,d]*10,self.cell[2,d]*10)
        s += "@AR3A\n"
        s += "{0}\n".format(len(self.reppositions))
        for i in range(len(self.reppositions)):
            position = np.dot(self.reppositions[i],self.cell)*10   #in Angstrom
            s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(position[0],
                                                    position[1],
                                                    position[2])
        s = "\n".join(self.doc) + "\n" + s
        print(s,end="")
