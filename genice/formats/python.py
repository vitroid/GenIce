from genice.formats.baseclass import GenIce


class Formatter(GenIce):
    """
    re-make python module for GenIce
    """
    def run(self, options):
        water_type    = options.water[0]
        guests        = options.guests
        self.stage1()   #replicate the unit cell
        self.stage2()   #prepare random graph
        self.stage3()   #Make an ice graph
        self.stage4()   #Depolarize
        self.stage5()   #Orientation
        self.stage6(water_type)  #Water atoms
        self.stage7(guests)      #Guest atoms
        self.logger.info("Total number of atoms: {0}".format(len(self.atoms)))
        self.logger.info("Output as a python module.")
        s = ""
        s += '"""\n'
        s += "\n".join(self.doc) + "\n"
        s += '"""\n'
        s += "bondlen={0}\n".format(self.bondlen)
        s += "coord='relative'\n"
        if self.cell[1,0] == 0 and self.cell[2,0] == 0 and self.cell[2,1] == 0:
            s += "celltype='rect'\n"
            s += "cell='{0} {1} {2}'\n".format(self.cell[0,0],self.cell[1,1],self.cell[2,2])
        else:
            s += "celltype='triclinic'\n"
            s += "cell='{0} {1} {2} {3} {4} {5} {6} {7} {8}'\n".format(*self.cell[0],*self.cell[1], *self.cell[2])
        s += "density={0}\n".format(self.density)
        s += "waters=\"\"\"\n"
        for i in range(len(self.reppositions)):
            position = self.reppositions[i]
            s += "{0:9.4f} {1:9.4f} {2:9.4f}\n".format(position[0],position[1],position[2])
        s += "\"\"\"\n\n"
        print(s,end="")
