import random
from collections import defaultdict
from logging import getLogger

import numpy as np
import pairlist as pl

from genice2 import digraph as dg
from genice2.cell import Cell
from genice2.decorators import banner, timeit
from genice2.genice import GenIce, put_in_array, shortest_distance
from genice2.valueparser import parse_pairs


class AnalIce(GenIce):
    def __init__(self, lat, signature=""):

        logger = getLogger()
        self.rep = (1, 1, 1)
        density = 0.0
        #         cations=dict(),
        #         anions=dict(),
        #         spot_guests=dict(),
        #         spot_groups=dict(),
        self.asis = False

        # Show the document of the module
        try:
            self.doc = lat.__doc__.splitlines()
        except BaseException:
            self.doc = []

        if len(signature) > 0:
            self.doc.append("")
            self.doc.append(signature)

        for line in self.doc:
            logger.info("  " + line)
        # ================================================================
        # rotmatrices (analice)
        #
        try:
            self.rotmatrices = lat.rotmat
        except BaseException:
            logger.info("No rotmatrices in lattice")
            self.rotmatrices = None
        # ================================================================
        # waters: positions of water molecules
        #
        self.waters = put_in_array(lat.waters)
        logger.debug("Waters: {0}".format(len(self.waters)))
        self.waters = self.waters.reshape((self.waters.size // 3, 3))

        # ================================================================
        # cell: cell dimension
        #   see parse_cell for syntax.
        #
        self.cell = Cell(lat.cell)

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative"
        #
        if lat.coord == "absolute":
            self.waters = self.cell.abs2rel(self.waters)

        self.waters = np.array([w - np.floor(w) for w in self.waters])

        # ================================================================
        # pairs: specify the pairs of molecules that are connected.
        #   Bond orientation will be shuffled later
        #   unless it is "fixed".
        #
        self.pairs = None

        try:
            self.pairs = parse_pairs(lat.pairs)
        except AttributeError:
            logger.info("HB connectivity is not defined.")

        # ================================================================
        # bondlen: specify the bond length threshold.
        #   This is used when "pairs" are not specified.
        #   It is applied to the original positions of molecules (before density setting).
        #
        nmol = self.waters.shape[0]  # nmol in a unit cell
        volume = self.cell.volume()  # volume of a unit cell in nm**3
        self.bondlen = None

        try:
            self.bondlen = lat.bondlen
            logger.info("Bond length (specified): {0}".format(self.bondlen))
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            # assume that the particles distribute homogeneously.
            rc = (volume / nmol)**(1 / 3) * 1.5
            p = pl.pairs_iter(self.waters,
                              rc=rc,
                              cell=self.cell.mat,
                              distance=False)
            self.bondlen = 1.1 * \
                shortest_distance(self.waters, self.cell, pairs=p)
            logger.info("Bond length (estim.): {0}".format(self.bondlen))

        # Set density
        mass = 18  # water
        NB = 6.022e23
        density0 = mass * nmol / (NB * volume * 1e-21)

        if density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                logger.info(
                    "Density is not specified. Assume the density from lattice.")
                dmin = shortest_distance(self.waters, self.cell)
                logger.info(
                    "Closest pair distance: {0} (should be around 0.276 nm)".format(dmin))
                self.density = density0 / (0.276 / dmin)**3
                # self.density = density0
        else:
            self.density = density

        logger.info("Target Density: {0}".format(self.density))
        logger.info("Original Density: {0}".format(density0))

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density)**(1.0 / 3.0)
        self.cell.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        logger.info("Bond length (scaled, nm): {0}".format(self.bondlen))

        # ================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        self.cagepos = None
        self.cagetype = None

        # ================================================================
        # fixed: specify the bonds whose directions are fixed.
        #   you can specify them in pairs at a time.
        #   You can also leave it undefined.
        #
        self.fixed = []
        try:
            self.fixed = parse_pairs(lat.fixed)
            logger.info("Orientations of some edges are fixed.")
        except AttributeError:
            pass

        self.dopeIonsToUnitCell = None
        self.dopants = set()

        # if asis, make pairs to be fixed.
        if self.asis and len(self.fixed) == 0:
            self.fixed = self.pairs

        # filled cages
        self.filled_cages = set()

        # groups info
        self.groups = defaultdict(dict)

    def analyze_ice(self, water, formatter, noise=0.):
        """
        Protocol for analice
        """

        def Stages():
            hooks = formatter.hooks()

            maxstage = max(0, *hooks.keys())

            if 0 in hooks:
                abort = hooks[0](self)
                if maxstage < 1 or abort:
                    return

            self.Stage1(noise)

            if 1 in hooks:
                abort = hooks[1](self)
                if maxstage < 2 or abort:
                    return

            if self.rotmatrices is None:
                res = self.Stage2()

            if 2 in hooks:
                abort = hooks[2](self)
                if maxstage < 3 or abort:
                    return

            if self.rotmatrices is None:
                self.Stage3()

            if 3 in hooks:
                abort = hooks[3](self)
                if maxstage < 4 or abort:
                    return

            self.Stage4()

            if 4 in hooks:
                abort = hooks[4](self)
                if maxstage < 5 or abort:
                    return

            # molecular orientation should be given in the loader.
            if self.rotmatrices is None:
                self.Stage5()

            if 5 in hooks:
                abort = hooks[5](self)
                if maxstage < 6 or abort:
                    return

            self.Stage6(water)

            if 6 in hooks:
                abort = hooks[6](self)
                if maxstage < 7 or abort:
                    return

            # self.Stage7_analice(guests)
            if 7 in hooks:
                hooks[7](self)

        abort = Stages()
        if not abort:
            return formatter.dump()

    @timeit
    @banner
    def Stage1(self,
               noise=0.):
        """Preparation.

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        """

        logger = getLogger()
        self.reppositions = self.waters

        # This must be done before the replication of the cell.
        logger.info("  Number of water molecules: {0}".format(
            len(self.reppositions)))

        # self.graph = self.prepare_random_graph(self.fixed)
        self.graph = self.prepare_random_graph(self.pairs)

        # scale the cell
        self.repcell = Cell(self.cell.mat)

        # self.repcell.scale2(self.rep)
        # add small perturbations to the molecular positions.
        if noise > 0.0:
            logger.info("  Add noise: {0}.".format(noise))
            perturb = np.random.normal(loc=0.0,
                                       scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                                       size=self.reppositions.shape)
            self.reppositions += self.repcell.abs2rel(perturb)

    @timeit
    @banner
    def Stage4(self):
        """Make a spacegraph.

        Provided variables:
        spacegraph: depolarized network with node positions.
        yapresult:  Animation of the depolarization process in YaPlot format.
        """

        logger = getLogger()
        self.yapresult = ""
        self.spacegraph = dg.SpaceIceGraph(self.graph,
                                           coord=self.reppositions,
                                           immutables=self.graph.immutables)
