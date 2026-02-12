# coding: utf-8

from logging import getLogger
import numpy as np
from genice3.genice import GenIce3
from io import TextIOWrapper
import sys

desc = {
    "ref": {},
    "brief": "Calculate the polarization of the ice.",
    "usage": "Calculate the polarization of the ice.",
}




def calculate(genice: GenIce3) -> np.ndarray:
    "Calculate the polarization of the ice."
    logger = getLogger("genice3.exporter._pol")

    polarization = np.zeros(3)
    for i,j in genice.digraph.edges():
        pos_i = genice.lattice_sites[i]
        pos_j = genice.lattice_sites[j]
        delta = pos_j - pos_i
        delta -= np.floor(delta + 0.5)
        polarization += delta
    return polarization

def dumps(genice: GenIce3, **options) -> str:
    """
    Dump the polarization of the ice.
    """
    polarization = calculate(genice)
    return f"{polarization[0]:.5f} {polarization[1]:.5f} {polarization[2]:.5f} polarization\n"


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **options):
    file.write(dumps(genice, **options))
