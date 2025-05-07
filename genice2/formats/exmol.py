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
        s += 'Lattice="{0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f} {8:.3f}"\n'.format(
            *genice.cell_matrix().reshape((9,)) * 10
        )
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:>4}{1:15.5f}{2:15.5f}{3:15.5f}\n".format(
                atomname, position[0] * 10, position[1] * 10, position[2] * 10
            )
        s = "#\n#Command line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
