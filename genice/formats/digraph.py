from genice.formats.baseclass import GenIce


class Formatter(GenIce):
    """
    Hydrogen bond network in @NGPH format.
    """
    def run(self, options):
        water_type    = options.water[0]
        guests        = options.guests
        self.stage1()   #replicate the unit cell
        self.stage2()   #prepare random graph
        self.stage3()   #Make an ice graph
        self.stage4()   #Depolarize
        self.logger.info("Output the hydrogen bond network.")
    
        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(self.reppositions))
        for i,j,k in self.graph.edges_iter(data=True):
            s += "{0} {1}\n".format(i,j)
        s += "-1 -1\n"
        s = "\n".join(self.doc) + "\n" + s
        print(s,end="")
