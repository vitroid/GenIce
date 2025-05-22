# coding: utf-8
"""
Crude Cif file format
"""

from logging import getLogger
from io import TextIOWrapper
import sys

import numpy as np

from genice2.decorators import timeit, banner
import genice2.formats
from genice2.molecules import serialize
from genice2.genice import GenIce
from genice2.cell import Cell


def nearly_zero(x):
    return x**2 < 1e-10


class Format(genice2.formats.Format):
    """
    Atomic positions are in a crude CIF format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    def format_cell_shape(self, genice: GenIce):
        aL, bL, cL, alpha, beta, gamma = Cell(genice.cell_matrix()).shape()

        s = ""
        s += "data_genice\n"
        s += "#" + "\n#Command line: " + " ".join(sys.argv) + "\n"
        s += "_cell_length_a                {0}\n".format(aL * 10)
        s += "_cell_length_b                {0}\n".format(bL * 10)
        s += "_cell_length_c                {0}\n".format(cL * 10)
        s += "_cell_angle_alpha             {0}\n".format(alpha)
        s += "_cell_angle_beta              {0}\n".format(beta)
        s += "_cell_angle_gamma             {0}\n".format(gamma)
        s += "\n"
        rights = np.array([alpha, beta, gamma]) - 90
        if np.allclose(rights, 0):
            s += "_symmetry_cell_setting        'orthorhombic'\n"
            s += "_symmetry_space_group_name_H-M   'P 1 '\n"
        else:
            # for now it is always triclinic
            s += "_symmetry_cell_setting        'triclinic'\n"
            s += "_symmetry_space_group_name_H-M   'P 1 '\n"
        return s

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output in CIF format."
        logger = getLogger()

        s = self.format_cell_shape(genice)

        s += """
loop_
  _symmetry_equiv_pos_as_xyz
X,Y,Z

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""
        universe = genice.full_atomic_positions()
        atoms = []
        for pos in universe:
            atoms += serialize(pos)

        cell = Cell(genice.cell_matrix())
        for i, atom in enumerate(atoms):
            molorder, resname, atomname, position, order = atom
            pos = cell.abs2rel(position)
            label = f"{atomname}{i}"
            s += f"{label:>6}{atomname:>6}{pos[0]:10.4f}{pos[1]:10.4f}{pos[2]:10.4f}\n"
        s += "\n"
        file.write(s)
