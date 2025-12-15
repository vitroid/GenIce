# coding: utf-8

from genice2.decorators import timeit, banner
import genice2.formats
from logging import getLogger
import numpy as np
desc = {
    "ref": {
        "K1939": 'J. G. Kirkwood, J. Chem. Phys. 7, 911 (1939).'},
    "brief": "Kirkworrd G factor.",
    "usage": "Kirkwood G is a convenient index for testing the long-range order.",
}


# It should be expressed as a function of distance.

class Format(genice2.formats.Format):
    """
Calculate the Kirkwood G factor.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def hooks(self):
        return {5: self.Hook5}

    @timeit
    @banner
    def Hook5(self, ice):
        "Kirkwood G."
        logger = getLogger()
        zvec = np.array([0., 0., 1.])
        dipoles = zvec @ ice.rotmatrices
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
        coord = ice.reppositions
        delta = np.zeros((n, n, 3))
        for i in (0, 1, 2):
            x = coord[:, i]
            xT = x.reshape((n, 1))
            d = x - xT
            d -= np.floor(d + 0.5)  # wrap
            delta[:, :, i] = d
        delta = delta @ ice.repcell.mat
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
        self.output = s
