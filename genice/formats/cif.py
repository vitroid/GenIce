from genice.formats.baseclass import GenIce
from math import acos, pi
import numpy as np

def nearly_zero(x):
    return np.dot(x,x) < 1e-10


class Formatter(GenIce):
    """
    Crude Cif file format
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
        self.logger.info("Output in MDView format.")

        a = self.cell[0,:]
        b = self.cell[1,:]
        c = self.cell[2,:]
        aL= np.linalg.norm(a)
        bL= np.linalg.norm(b)
        cL= np.linalg.norm(c)
        ab = np.dot(a,b)
        bc = np.dot(b,c)
        ca = np.dot(c,a)
        alpha = acos(bc/(bL*cL)) * 180 / pi
        beta  = acos(ca/(cL*aL)) * 180 / pi
        gamma = acos(ab/(aL*bL)) * 180 / pi
        s = ""
        #if celltype == "rect":
        #    s += "-length '({0}, {1}, {2})'\n".format(cell[0,0]*10,cell[1,1]*10,cell[2,2]*10)
        s += "data_geniice_{0}\n".format(options.Type[0])
        s += '#' + "\n#".join(self.doc) + "\n"
        s += "_cell_length_a                {0}\n".format(aL*10)
        s += "_cell_length_b                {0}\n".format(bL*10)
        s += "_cell_length_c                {0}\n".format(cL*10)
        s += "_cell_angle_alpha             {0}\n".format(alpha)
        s += "_cell_angle_beta              {0}\n".format(beta)
        s += "_cell_angle_gamma             {0}\n".format(gamma)
        s += "\n"
        if nearly_zero(alpha-90) and nearly_zero(beta-90) and nearly_zero(gamma-90):
            s += "_symmetry_cell_setting        'orthorhombic'"
        else:
            s += "_symmetry_cell_setting        'triclinic'"  #for now it is always triclinic
        s += """
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""
        celli = np.linalg.inv(self.cell)
        for i in range(len(self.atoms)):
            molorder, resname, atomname, position = self.atoms[i]
            position = np.dot(position, celli)
            label = atomname+"{0}".format(i)
            s += "{4:>6}{0:>6}{1:10.4f}{2:10.4f}{3:10.4f}\n".format(atomname,position[0],position[1],position[2],label)
        s += "\n"
        print(s,end="")
