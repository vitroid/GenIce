# coding: utf-8
"""
A formatter plugin for GenIce2 to write in LAMMPS format.

Usage:
    genice2 ice5 -f lammps[mW] > ice5.lammps

Options:
    mW   mW water model. (No guests are available).

Note: all-atom model is experimental.
"""

from math import cos, radians
from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from logging import getLogger
import numpy as np
import sys

desc = {"ref": {},
        "brief": "LAMMPS file (TESTING).",
        "usage": __doc__
        }


AA = 0.1  # nm


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in LAMMPS format.
No options available.
    """

    def __init__(self, **kwargs):
        self.mW = False
        unknown = dict()
        for k, v in kwargs.items():
            if k == "mW" and v is True:
                # for commandline use
                self.mW = True
            else:
                unknown[k] = v
        super().__init__(**kwargs)

    def hooks(self):
        return {1:self.Hook1, 7: self.Hook7}

    def _header(self, natoms):
        s = "# Generated by GenIce2\n"
        s += "\n" # this newline is necessary
        s += f"  {natoms} atoms\n"
        return s

    def _cellshape(self, cellmat):
        a,b,c,A,B,C = genice2.cell.cellshape(cellmat)
        # https://docs.lammps.org/Howto_triclinic.html
        # I hope it works.
        lx = a
        xy = b*cos(radians(C))
        xz = c*cos(radians(B))
        ly = (b**2-xy**2)**0.5
        yz = (b*c*cos(radians(A)) - xy*xz)/ly
        lz = (c**2 - xz**2 - yz**2)**0.5
        # logger.info([a,b,c,A,B,C])
        # logger.info(cellmat)
        # logger.info([lx,ly,lz,xy,xz,yz])

        s = ""
        s += f"0.0 {lx} xlo xhi\n"
        s += f"0.0 {ly} ylo yhi\n"
        s += f"0.0 {lz} zlo zhi\n"
        if not -1e10 < xy < 1e-10:
            s += f"{xy} xy\n"
        if not -1e10 < xz < 1e-10:
            s += f"{xz} xz\n"
        if not -1e10 < yz < 1e-10:
            s += f"{yz} yz\n"
        s += "\n"
        return s

    @timeit
    @banner
    def Hook1(self, ice):
        "Output in LAMMPS format."
        if not self.mW:
            return
        logger = getLogger()

        conv = 1.0 / AA

        atoms = ice.reppositions @ ice.repcell.mat * conv
        logger.info("  Total number of atoms: {0}".format(len(atoms)))

        s = self._header(len(atoms))
        s += f"  1 atom types\n"
        s += "\n"

        s += self._cellshape(ice.repcell.mat / AA)

        s += "Masses\n"
        s += "\n"
        masses = {"mW":18.015, }
        for i, atom in enumerate(masses):
            s += f" {i+1} {masses[atom]} # {atom}\n"
        s += "\n"
        s += "Atoms\n"
        s += "\n"
        for i, atom in enumerate(atoms):
            s += f"{i+1} {i+1} 1 0 {atom[0]} {atom[1]} {atom[2]}\n"
        self.output = s
        return True

    @timeit
    @banner
    def Hook7(self, ice):
        "Output in LAMMPS format."
        logger = getLogger()
        logger.warning("All-atom model is experimental.")

        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)
        logger.info("  Total number of atoms: {0}".format(len(atoms)))
        s = self._header(len(atoms))
        # atom types
        atomdic = dict()
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            if atomname not in atomdic:
                atomdic[atomname] = len(atomdic) + 1
        s += f"  {len(atomdic)} atom types\n"
        s += "\n"

        s += self._cellshape(ice.repcell.mat / AA)

        s += "Masses\n"
        s += "\n"
        masses = {"O":15.999, "H":1.0, "C":12.0, "M":0.0}
        for atom in sorted(atomdic, key=lambda x:atomdic[x]):
            s += f" {atomdic[atom]} {masses[atom[0]]} # {atom}\n"
        s += "\n"
        s += "Atoms\n"
        s += "\n"
        for i, atom in enumerate(atoms):
            molorder, resname, atomname, position, order = atom
            s += f"{i+1} {atomdic[atomname]} {position[0] / AA} {position[1] / AA} {position[2] / AA}\n"
        self.output = s
