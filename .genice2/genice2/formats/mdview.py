# coding: utf-8

from logging import getLogger
from io import TextIOWrapper

from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from genice2.genice import GenIce

desc = {
    "ref": {},
    "brief": "MDView file (in Angdtrom).",
    "usage": "No options available.",
    "test": (
        {"args": "au"},
        {"args": "AA"},
        {"args": {"au": True}},
        {"args": {"AA": True}},
    ),
}

if __name__[-6:] == "mdv_au":
    desc["brief"] = "MDView file (in au)."


au = 0.052917721092  # nm
AA = 0.1  # nm


class Format(genice2.formats.Format):
    """
    The atomic positions of the molecules are output in MDView format.
    No options available.
    """

    def __init__(self, **kwargs):
        if __name__[-6:] == "mdv_au":
            self.conv = 1.0 / au
        else:
            self.conv = 1.0 / AA
        for k, v in kwargs.items():
            assert v == True and k in ("au", "AA"), f"Unknown keyword: {k}"
            if k == "au":
                self.conv = 1.0 / au
            else:
                self.conv = 1.0 / AA

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output in MDView format."
        logger = getLogger()
        atoms = []
        for mols in genice.full_atomic_positions():
            atoms += serialize(mols)
        logger.info("  Total number of atoms: {0}".format(len(atoms)))
        cellmat = genice.cell_matrix()
        s = ""
        s += "# {0} {1} {2}\n".format(
            cellmat[0, 0] * self.conv,
            cellmat[1, 1] * self.conv,
            cellmat[2, 2] * self.conv,
        )
        s += "-center 0 0 0\n"
        s += "-fold\n"
        s += "{0}\n".format(len(atoms))
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            s += "{0:5} {1:9.4f} {2:9.4f} {3:9.4f}\n".format(
                atomname, *(position[:3] * self.conv)
            )
        file.write(s)
