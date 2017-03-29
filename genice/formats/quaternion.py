from genice.formats.baseclass import GenIce
from genice import rigid
import numpy     as np


class Formatter(GenIce):
    """
    Rigid water molecule 
    """
    def run(self, options):
        water_type    = options.water[0]
        guests        = options.guests
        self.stage1()   #replicate the unit cell
        self.stage2()   #prepare random graph
        self.stage3()   #Make an ice graph
        self.stage4()   #Depolarize
        self.stage5()   #Orientation
        self.logger.info("Output water molecules as rigid rotors (Quaternion).")
        s = ""
        if self.cell[1,0] == 0 and self.cell[2,0] == 0 and self.cell[2,1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(self.cell[0,0]*10,self.cell[1,1]*10,self.cell[2,2]*10)
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(self.cell[0,d]*10,self.cell[1,d]*10,self.cell[2,d]*10)
        s += "@NX4A\n"
        s += "{0}\n".format(len(self.reppositions))
        for i in range(len(self.reppositions)):
            position = np.dot(self.reppositions[i],self.cell)*10   #in Angstrom
            quat     = rigid.rotmat2quat(self.rotmatrices[i].transpose())
            s += "{0:9.4f} {1:9.4f} {2:9.4f}  {3:9.4f} {4:9.4f} {5:9.4f} {6:9.4f}\n".format(position[0],
                                                                                position[1],
                                                                                position[2],
                                                                                quat[0],
                                                                                quat[1],
                                                                                quat[2],
                                                                                quat[3])
        s = "\n".join(self.doc) + "\n" + s
        print(s,end="")
