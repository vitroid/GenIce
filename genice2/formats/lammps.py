# coding: utf-8

from genice2.molecules import serialize
import genice2.formats
from genice2.decorators import timeit, banner
from logging import getLogger
import numpy as np
import sys

desc = {"ref": {},
        "brief": "LAMMPS file (in Angdtrom).",
        "usage": "No options available."
        }


AA = 0.1  # nm


class Format(genice2.formats.Format):
    """
The atomic positions of the molecules are output in LAMMPS format.
No options available.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {7: self.Hook7}

    @timeit
    @banner
    def Hook7(self, ice):
        "Output in LAMMPS format."
        logger = getLogger()
        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)
        logger.info("  Total number of atoms: {0}".format(len(atoms)))
        conv = 1.0 / AA
        cellmat = ice.repcell.mat
        c2 = cellmat.copy()
        c2[0,0] = 0.0
        c2[1,1] = 0.0
        c2[2,2] = 0.0
        if not np.allclose(c2, 0.0):
            logger.error("LAMMPS plugin does not accept a non-cuboidal cell.")
            sys.exit(1)
        s = "# Generated by GenIce2\n"
        s += "\n" # this newline is necessary
        s += f"  {len(atoms)} atoms\n"

        # atom types
        atomdic = dict()
        for atom in atoms:
            molorder, resname, atomname, position, order = atom
            if atomname not in atomdic:
                atomdic[atomname] = len(atomdic) + 1
        s += f"  {len(atomdic)} atom types\n"
        s += "\n"
        s += f"0.0 {cellmat[0,0]*conv} xlo xhi\n"
        s += f"0.0 {cellmat[1,1]*conv} ylo yhi\n"
        s += f"0.0 {cellmat[2,2]*conv} zlo zhi\n"
        s += "\n"
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
            s += f"{i+1} {atomdic[atomname]} {position[0]*conv} {position[1]*conv} {position[2]*conv}\n"
        self.output = s
