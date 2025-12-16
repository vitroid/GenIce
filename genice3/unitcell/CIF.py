# It is a dummy unitcell to load a cif file.
# usage: genice3 CIF[file=diamond.cif, site=C]


import genice3.unitcell
import numpy as np
import networkx as nx
from cif2ice import cellvectors, read_cif
from genice3.util import atoms_to_waters, shortest_distance
import re
from logging import getLogger


class UnitCell(genice3.unitcell.UnitCell):
    """
    cif単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        logger = getLogger("CIF")
        file = kwargs.get("file")
        osite = kwargs.get("osite")
        hsite = kwargs.get("hsite")
        if file is None:
            raise ValueError("file is required")
        if osite is None:
            osite = "O"

        # download(URL, fNameIn)

        atoms, box = read_cif.read_and_process(file, make_rect_box=False)
        # pattern matching
        oatoms = np.array([a[1:] for a in atoms if re.match(osite, a[0])])
        logger.info(f"{osite=} {oatoms.shape=} {atoms}")
        cell = cellvectors(*box) / 10  # nm
        shortest_OO = shortest_distance(oatoms, cell)
        cell *= 0.276 / shortest_OO

        if hsite is None:
            # 水素位置は指定されていないので、genice3にまかせる。
            super().__init__(
                cell=cell,  # nm
                waters=oatoms,
                coord="relative",
            )
        else:
            hatoms = np.array([a[1:] for a in atoms if re.match(hsite, a[0])])
            waters, pairs, oo_pairs = atoms_to_waters(
                oatoms, hatoms, cell, partial_order=True
            )
            super().__init__(
                cell=cell,
                waters=waters,
                graph=nx.Graph(oo_pairs),
                fixed=nx.DiGraph(pairs),
                coord="relative",
            )
