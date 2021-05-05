# coding: utf-8
import numpy as np
from logging import getLogger
import genice2.molecules

desc = {
    "ref": {"NvdE": "Nada 2003"},
    "usage": "No options available.",
    "brief": "A 6-site water model.",
}
water = 1  # Identify


class Molecule(genice2.molecules.Molecule):

    def __init__(self):
        lx = 0.732813006922696 / 10
        hy = 0.792836654487448 / 10
        oz = -6.400328302740263E-002 / 10
        hz = 0.512026264219221 / 10
        mz = 0.165996716972597 / 10
        lz = -0.567651708900964 / 10
        self.sites_ = np.array([[0, 0, oz],
                                [0, hy, hz],
                                [0, -hy, hz],
                                [0, 0, mz],
                                [lx, 0, lz],
                                [-lx, 0, lz],
                                ])  # nm, OHHMLL

        self.labels_ = ["O", "H", "H", "M", "L", "L"]
        self.name_ = "SOL"


#  0.000000000000000E+000  0.000000000000000E+000 -6.400328302740263E-002 16 O
#  0.000000000000000E+000  0.792836654487448       0.512026264219221      1  H
#  0.000000000000000E+000 -0.792836654487448       0.512026264219221      1  H
#  0.000000000000000E+000  0.000000000000000E+000  0.165996716972597      0  M
#  0.732813006922696       0.000000000000000E+000 -0.567651708900964      0  L
# -0.732813006922696       0.000000000000000E+000 -0.567651708900964      0  L
