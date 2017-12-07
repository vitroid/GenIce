# coding: utf-8
"""
Crude Cif file format
"""

from math import acos,pi
import numpy as np

def nearly_zero(x):
    return np.dot(x,x) < 1e-10


def hook7(lattice):
    lattice.logger.info("Hook7: Output in CIF format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    cellmat = lattice.repcell.mat
    
    a = cellmat[0,:]
    b = cellmat[1,:]
    c = cellmat[2,:]
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
    s += "data_genice_{0}\n".format(lattice.lattice_type)
    s += '#' + "\n#".join(lattice.doc) + "\n"
    s += "_cell_length_a                {0}\n".format(aL*10)
    s += "_cell_length_b                {0}\n".format(bL*10)
    s += "_cell_length_c                {0}\n".format(cL*10)
    s += "_cell_angle_alpha             {0}\n".format(alpha)
    s += "_cell_angle_beta              {0}\n".format(beta)
    s += "_cell_angle_gamma             {0}\n".format(gamma)
    s += "\n"
    if nearly_zero(alpha-90) and nearly_zero(beta-90) and nearly_zero(gamma-90):
        s += "_symmetry_cell_setting        'orthorhombic'\n"
        s += "_symmetry_space_group_name_H-M   'P 1 '\n"
    else:
        s += "_symmetry_cell_setting        'triclinic'\n"  #for now it is always triclinic
        s += "_symmetry_space_group_name_H-M   'P 1 '\n"
    s += """
loop_
  _symmetry_equiv_pos_as_xyz
X,Y,Z

loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
"""
    for i in range(len(lattice.atoms)):
        molorder, resname, atomname, position, order = lattice.atoms[i]
        position = lattice.repcell.abs2rel(position)
        label = atomname+"{0}".format(i)
        s += "{4:>6}{1:10.4f}{2:10.4f}{3:10.4f}{0:>6}\n".format(atomname,position[0],position[1],position[2],label)
    s += "\n"
    print(s,end="")
    lattice.logger.info("Hook7: end.")

hooks = {7:hook7}
