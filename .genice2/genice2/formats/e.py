# coding: utf-8

from logging import getLogger
from io import TextIOWrapper
import sys

from genice2 import rigid
import genice2.formats
from genice2.decorators import timeit, banner
from genice2.genice import GenIce

desc = {
    "ref": {},
    "brief": "Rigid rotor (Euler angle).",
    "usage": "No options.",
}


class Format(genice2.formats.Format):
    """
    The positions and orientations (in Euler angles) of the water molecules are output in @NX3A format.
    No options available.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Output water molecules as rigid rotors (Euler)."
        logger = getLogger()
        cellmat = genice.cell_matrix()
        s = ""
        if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
            s += "@BOX3\n"
            s += "{0} {1} {2}\n".format(
                cellmat[0, 0] * 10, cellmat[1, 1] * 10, cellmat[2, 2] * 10
            )
        else:
            s += "@BOX9\n"
            for d in range(3):
                s += "{0} {1} {2}\n".format(
                    cellmat[0, d] * 10, cellmat[1, d] * 10, cellmat[2, d] * 10
                )
        water_positions = genice.water_positions()
        rotation_matrices = genice.rotation_matrices()

        s += "@NX3A\n"
        s += "{0}\n".format(len(water_positions))
        for pos, rot in zip(water_positions, rotation_matrices):
            position = pos @ cellmat * 10  # in Angstrom
            euler = rigid.quat2euler(rigid.rotmat2quat(rot.transpose()))
            s += "{0:9.4f} {1:9.4f} {2:9.4f}  {3:9.4f} {4:9.4f} {5:9.4f}\n".format(
                position[0], position[1], position[2], euler[0], euler[1], euler[2]
            )
        s = "\n# Command line: " + " ".join(sys.argv) + "\n" + s
        file.write(s)
