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


@dataclass(frozen=True)
class AtomicStructure:
    resname: str
    atomnames: list[str]
    positions: list[np.ndarray]
    orig_order: list[int]

    @classmethod
    def from_molecule(
        cls,
        rel_coords: np.ndarray,
        cell: Cell,
        rotmatrices: np.ndarray,
        molecule: Molecule,
        immutables=set(),
    ):
        logger = getLogger()
        logger.info(f"from_molecule: {molecule}")
        cls.resname = molecule.name
        cls.atomnames = molecule.labels
        cls.positions = []
        cls.orig_order = []

        for order, rel_coord in enumerate(rel_coords):
            if order in immutables:
                continue
            cls.orig_order.append(order)

            abscom = cell.rel2abs(rel_coord)  # relative to absolute
            rotated = molecule.sites @ rotmatrices[order]
            cls.positions.append(abscom + rotated)

        return cls

    @classmethod
    def monatom(cls, rel_coord: np.ndarray, cell: Cell, name: str):
        logger = getLogger()

        cls.resname = name
        cls.atomnames = [name]
        cls.positions = [
            np.array(
                [
                    cell.rel2abs(rel_coord),
                ]
            )
        ]
        cls.orig_order = [0]

        return cls


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
