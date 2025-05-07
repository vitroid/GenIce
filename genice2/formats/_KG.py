# coding: utf-8

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
import numpy as np
from io import TextIOWrapper
from genice2.genice import GenIce

desc = {
    "ref": {"K1939": "J. G. Kirkwood, J. Chem. Phys. 7, 911 (1939)."},
    "brief": "Kirkworrd G factor.",
    "usage": "Kirkwood G is a convenient index for testing the long-range order.",
}


# It should be expressed as a function of distance.


class Format(genice2.formats.Format):
    """
    Calculate the Kirkwood G factor.
    """

    def __init__(self, **kwargs):
        pass

    @timeit
    @banner
    def dump(self, genice: GenIce, file: TextIOWrapper):
        "Kirkwood G."
        logger = getLogger()
        zvec = np.array([0.0, 0.0, 1.0])
        dipoles = zvec @ genice.rotation_matrices()
        n = dipoles.shape[0]
        logger.info("  Inner products of the dipoles.")
        # inner products of dipole pairs
        ip = np.zeros((n, n, 3))
        for i in (0, 1, 2):  # x,y,z
            x = dipoles[:, i]
            xT = x.reshape((n, 1))
            ip[:, :, i] = x * xT  # broadcast product
        ip = np.sum(ip, axis=2).reshape(n * n)
        logger.info("  Pair distances.")
        # pair distances
        coord = genice.water_positions()
        delta = np.zeros((n, n, 3))
        for i in (0, 1, 2):
            x = coord[:, i]
            xT = x.reshape((n, 1))
            d = x - xT
            d -= np.floor(d + 0.5)  # wrap
            delta[:, :, i] = d
        delta = delta @ genice.cell_matrix()
        delta = delta * delta
        delta = np.sum(delta, axis=2)
        delta = np.sqrt(delta).reshape(n * n)
        # print(ip)
        # print(delta)
        # print(ice.repcell.mat)
        logger.info("  G(r)")
        mx = np.max(delta)
        x = 0.04
        s = ""
        while x < mx:
            filter = ip[delta < x]
            s += "{0:.5f} {1:.5f}\n".format(x, np.sum(filter) / n)
            x += 0.08
        file.write(s)
