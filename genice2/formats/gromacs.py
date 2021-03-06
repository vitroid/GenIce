# coding: utf-8

desc={"ref": {"gro": "http://manual.gromacs.org/current/online/gro.html"},
      "brief": "Gromacs .gro file.",
      "usage": "No options available."
      }


from logging import getLogger
import genice2.formats
from genice2.decorators import timeit, banner


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in Gromacs format.
No options available.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {7:self.Hook7}


    @timeit
    @banner
    def Hook7(self, ice):
        "Output in Gromacs format."
        logger = getLogger()
        logger.info("  Total number of atoms: {0}".format(len(ice.atoms)))
        if len(ice.atoms) > 99999:
            logger.warn("  Gromacs fixed format cannot deal with atoms more than 99999. Residue number and atom number are faked.")
        cellmat = ice.repcell.mat
        s = ""
        s += "Generated by GenIce https://github.com/vitroid/GenIce \n"
        s += "{0}\n".format(len(ice.atoms))
        molorder = 0
        for i, atom in enumerate(ice.atoms):
            resno, resname, atomname, position, order = atom
            if resno == 0:
                molorder += 1
            if len(ice.atoms) > 99999:
                s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(9999,resname, atomname, 9999,position[0],position[1],position[2])
            else:
                s += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(molorder,resname, atomname, i+1,position[0],position[1],position[2])
        if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
            s += "    {0} {1} {2}\n".format(cellmat[0,0],cellmat[1,1],cellmat[2,2])
        else:
            assert cellmat[0,1] == 0 and cellmat[0,2] == 0 and cellmat[1,2] == 0
            s += "    {0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(cellmat[0,0],
                                                                    cellmat[1,1],
                                                                    cellmat[2,2],
                                                                    cellmat[0,1],
                                                                    cellmat[0,2],
                                                                    cellmat[1,0],
                                                                    cellmat[1,2],
                                                                    cellmat[2,0],
                                                                    cellmat[2,1],
                                                                    )
        s += '#' + "\n#".join(ice.doc) + "\n"
        self.output = s
