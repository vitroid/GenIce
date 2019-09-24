# coding: utf-8

desc={"ref": {"II": 'Nakamura, Tatsuya et al. “Thermodynamic Stability of Ice II and Its Hydrogen-Disordered Counterpart: Role of Zero-Point Energy.” The Journal of Physical Chemistry B 120.8 (2015): 1843–1848.'},
         "brief": "Generate a variant of ice II by shuffling homodromic rings.",
         "usage": "-f nshuffle[x]\nx\tFraction of inverted hexagonal rings.",
         }

from genice.formats.euler import hook5
from countrings import countrings_nx as cr
from logging import getLogger
from attrdict import AttrDict
import re
import random

options = AttrDict()

def hook3(lattice):
    """
    It is called after defects are eliminated from the graph.
    """
    logger = getLogger()
    logger.info("Hook 3: Invert some hexagonal rings.")
    logger = getLogger()
    rings = [ring for ring in cr.CountRings(lattice.graph).rings_iter(6)
             if len(ring) == 6 and
             lattice.graph.is_homodromic(ring, cyclic=True)]
    Nr = len(rings)
    logger.info("  Nring:{0}".format(Nr))
    Ninv = Nr*options.frac//100
    invs = set()
    while len(invs) < Ninv:
        invs.add(random.randint(0,Nr-1))
    logger.info("  Invert: {0}".format(invs))
    for ringid in invs:
        lattice.graph.invert_path(rings[ringid], cyclic=True, force=True)
    logger.info("Hook 3: end.")
    
def argparser(self, arg):
    logger = getLogger()
    logger.info("Hook0: Options")
    assert re.match("^[0-9]+$", arg) is not None, "Argument must be an integer."
    options.frac = int(arg)
    logger.info("  Fraction of homodromic cycles to be inverted: {0}".format(options.frac))
    logger.info("Hook0: end.")

hooks = {0:argparser, 3:hook3, 5:hook5}
