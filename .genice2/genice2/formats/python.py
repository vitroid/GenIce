# coding: utf-8

import itertools as it
import re
from logging import getLogger
from math import floor

import genice2.formats
import numpy as np
from genice2.decorators import banner, timeit

desc = {
    "ref": {},
    "brief": "Cell-reshaper.",
    "usage": """
A formatter plugin for GenIce to produce a python lattice plugin.

Usage:
    genice2 ice5 --reshape 3,0,0,1,1,0,0,0,1 -f python > reshaped_ice5.py

Options:
    None.
""",
}


class Format(genice2.formats.Format):
    """
    A formatter plugin to produce a python lattice plugin.

    Options:
        None.
    """

    def __init__(self, **kwargs):
        logger = getLogger()
        self.ijk = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        unknown = dict()
        for k, v in kwargs.items():
            unknown[k] = v
        if len(unknown) > 0:
            logger.warning("Options in the python-formatting plugin is deprecated.")
        super().__init__(**unknown)

    def hooks(self):
        return {1: self.Hook1}

    @timeit
    @banner
    def Hook1(self, ice):
        "Make a python module."
        # Original cell matrix.
        cellmat = ice.repcell.mat

        # header
        s = ""
        s += '"""\n'
        s += "\n".join(ice.doc) + "\n"
        s += '"""\n'
        s += "import genice2.lattices\n\n"

        s += "class Lattice(genice2.lattices.Lattice):\n"
        s += "    def __init__(self):\n"
        s += "        self.bondlen={0}\n".format(ice.bondlen)
        s += "        self.coord='relative'\n"
        if (
            np.allclose(cellmat[1, 0], 0)
            and np.allclose(cellmat[2, 0], 0)
            and np.allclose(cellmat[2, 1], 0)
        ):
            s += "        from genice2.cell import cellvectors\n"
            s += f"        self.cell = cellvectors(a={abs(cellmat[0, 0]):.8f}, b={abs(cellmat[1, 1]):.8f}, c={abs(cellmat[2, 2]):.8f})\n"
        else:
            s += "\n        import numpy as np\n"
            s += "        self.cell=np.array(["
            for d in range(3):
                s += f"        [{cellmat[d, 0]:.8f}, {cellmat[d, 1]:.8f}, {cellmat[d, 2]:.8f}], "
            s += "        ])\n"
        s += "        self.density={0}\n".format(ice.density)
        s += "        self.waters="

        ss = '"""' + "\n".join(
            [" ".join([f"{y:.4f}" for y in x]) for x in ice.reppositions]
        )

        s += ss + '\n"""' + "\n"

        # 一時的に不活性にする。
        # # if ice.cagepos1 is not None:
        #     s += '        self.cages="""' + "\n"
        #     ncell, ss = FindEmptyCells(
        #         cellmat, self.ijk, ice.repcagepos, labels=ice.repcagetype
        #     )
        #     s += ss + '"""' + "\n\n"

        self.output = s
