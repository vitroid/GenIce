# coding: utf-8
from logging import getLogger
from io import TextIOWrapper
import sys

from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from genice2.genice import GenIce

desc = {
    "ref": {
        "exyz": "https://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html"
    },
    "brief": "Extended XYZ format.",
    "usage": "No options available.",
}


class Format(genice2.formats.Format):
    """
    The atomic positions of the molecules are output in an extended XYZ format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output in extended XYZ format."
        logger = getLogger()
        atoms = []
        for mols in genice.full_atomic_positions():
            atoms += serialize(mols)
        s = ""
        s += f"{len(atoms)}\n"
        s += "%PBC\n"
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            s += f"{atomname:>4}{position[0] * 10:15.5f}{position[1] * 10:15.5f}{position[2] * 10:15.5f}\n"
        s = "#" + "\n#Command line: " + " ".join(sys.argv) + "\n" + s
        s += "\n"
        cell = genice.cell_matrix() * 10
        s += f"Vector1 {cell[0][0]:15.5f}{cell[0][1]:15.5f}{cell[0][2]:15.5f}\n"
        s += f"Vector2 {cell[1][0]:15.5f}{cell[1][1]:15.5f}{cell[1][2]:15.5f}\n"
        s += f"Vector3 {cell[2][0]:15.5f}{cell[2][1]:15.5f}{cell[2][2]:15.5f}\n"
        s += f"Offset  {0.0:15.5f}{0.0:15.5f}{0.0:15.5f}\n"
        file.write(s)
