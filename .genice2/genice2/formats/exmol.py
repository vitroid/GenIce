# coding: utf-8
from logging import getLogger
from io import TextIOWrapper
import sys

from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from genice2.genice import GenIce

desc = {
    "ref": {"exmol": "http://libatoms.github.io/QUIP/io.html#extendedxyz"},
    "brief": "Extended XMol file format.",
    "usage": "No options available.",
}


class Format(genice2.formats.Format):
    """
    The atomic positions of the molecules are output in XMol format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output in extended XMol file format."
        logger = getLogger()
        atoms = []
        for mols in genice.full_atomic_positions():
            atoms += serialize(mols)
        s = ""
        s += "{0}\n".format(len(atoms))
        cell = genice.cell_matrix().reshape((9,)) * 10
        s += f'Lattice="{cell[0]:.3f} {cell[1]:.3f} {cell[ 2]:.3f} {cell[3]:.3f} {cell[4]:.3f} {cell[5]:.3f} {cell[6]:.3f} {cell[7]:.3f} {cell[8]:.3f}"\n'
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            s += f"{atomname:>4}{position[0] * 10:15.5f}{position[1] * 10:15.5f}{position[2] * 10:15.5f}\n"
        s = "#\n#Command line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
