# coding: utf-8
from logging import getLogger
from dataclasses import dataclass, field
import numpy as np

from genice2.cell import Cell


@dataclass
class Molecule:
    """
    Base class of a molecule
    """

    sites: np.ndarray
    labels: list[str]
    name: str
    is_water: bool = False


class AtomicStructure:
    def __init__(
        self,
        rel_coords: np.ndarray,
        cell: Cell,
        rotmatrices: np.ndarray = None,
        molecule: Molecule = None,
        immutables: set = None,
        name: str = None,
    ):
        logger = getLogger()

        if rotmatrices is None:
            # a shortcut for monatomic molecules
            self.resname = name
            self.atomnames = [name] * len(rel_coords)
            self.positions = cell.rel2abs(rel_coords).reshape(-1, 1, 3)
            self.orig_order = [0] * len(rel_coords)
            return

        self.resname = molecule.name
        self.atomnames = molecule.labels
        self.positions = []
        self.orig_order = []

        for order, rel_coord in enumerate(rel_coords):
            if order in immutables:
                continue
            self.orig_order.append(order)

            abscom = cell.rel2abs(rel_coord)  # relative to absolute
            rotated = molecule.sites @ rotmatrices[order]
            self.positions.append(abscom + rotated)


def serialize(mols: AtomicStructure):
    """
    mols is a collection of molecules
    make them into a list of atoms for Gromacs
    """
    logger = getLogger()
    atoms = []
    for i, atompositions in enumerate(mols.positions):
        for j, atompos in enumerate(atompositions):
            atoms.append(
                [j, mols.resname, mols.atomnames[j], atompos, mols.orig_order[i]]
            )
    return atoms
