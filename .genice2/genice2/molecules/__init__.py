# coding: utf-8
from logging import getLogger
from types import SimpleNamespace

import numpy as np


def arrange(coord, cell, rotmatrices, molecule, immutables=set()):
    logger = getLogger()

    name, labels, intra = molecule.get()  # Molecule class
    mols = SimpleNamespace(resname=name,
                           atomnames=labels,
                           positions=[],         # atomic positions
                           orig_order=[],  #
                           )
    for order, pos in enumerate(coord):
        if order in immutables:
            continue
        mols.orig_order.append(order)

        abscom = cell.rel2abs(pos)  # relative to absolute
        rotated = intra @ rotmatrices[order]
        mols.positions.append(abscom + rotated)
    return mols


def monatom(coord, cell, name):
    logger = getLogger()

    mols = SimpleNamespace(resname=name,
                           atomnames=[name],
                           # atomic positions
                           positions=[np.array([cell.rel2abs(coord), ])],
                           orig_order=[0],  #
                           )
    return mols


def serialize(mols):
    """
    mols is a collection of molecules
    make them into a list of atoms for Gromacs
    """
    logger = getLogger()
    atoms = []
    for i, atompositions in enumerate(mols.positions):
        for j, atompos in enumerate(atompositions):
            atoms.append([j, mols.resname, mols.atomnames[j],
                         atompos, mols.orig_order[i]])
    return atoms


class Molecule():
    """
    Base class of molecules
    """

    def __init__(self, **kwargs):
        assert len(kwargs) == 0

        # sites: positions of interaction sites relative to a center of molecule
        # a numpy array of rank (N, 3) where N is number of sites
        self.sites_ = np.zeros([1, 3])

        # Labels of the interaction sites.
        self.labels_ = ["Me", ]

        # the name that represents the molecule. It is necessary for Gromacs
        # format.
        self.name_ = "MET"

    def get(self):
        """
        Return an instance of the molecule.
        """
        return self.name_, self.labels_, self.sites_
