import logging
import sys
import numpy as np
from genice.cell   import Cell
from genice.genice import GenIce, put_in_array, shortest_distance
from genice.valueparsers import parse_pairs
import pairlist as pl
from collections import defaultdict
from genice import digraph as dg

class AnalIce(GenIce):
    def __init__(self, lat):

        self.logger = logging.getLogger()
        self.rep = (1,1,1)
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

        self.doc.append("")
        self.doc.append("Command line: {0}".format(" ".join(sys.argv)))

        for line in self.doc:
            self.logger.info("  "+line)
        # ================================================================
        # rotmatrices (analice)
        #
        try:
            self.rotmatrices = lat.rotmat
        except BaseException:
            self.logger.info("No rotmatrices in lattice")
            self.rotmatrices = None
        # ================================================================
        # waters: positions of water molecules
        #
        self.waters = put_in_array(lat.waters)
        self.logger.debug("Waters: {0}".format(len(self.waters)))
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
            self.logger.info("HB connectivity is not defined.")

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
            self.logger.info("Bond length (specified): {0}".format(self.bondlen))
        except AttributeError:
            self.logger.debug("  Estimating the bond threshold length...")
            # assume that the particles distribute homogeneously.
            rc = (volume / nmol)**(1 / 3) * 1.5
            grid = pl.determine_grid(self.cell.mat, rc)
            p = pl.pairs_fine(self.waters, rc, self.cell.mat, grid, distance=False)
            self.bondlen = 1.1 * shortest_distance(self.waters, self.cell, pairs=p)
            self.logger.info("Bond length (estim.): {0}".format(self.bondlen))

        # Set density
        mass = 18  # water
        NB = 6.022e23
        density0 = mass * nmol / (NB * volume * 1e-21)

        if density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                self.logger.info(
                    "Density is not specified. Assume the density from lattice.")
                dmin = shortest_distance(self.waters, self.cell)
                self.logger.info(
                    "Closest pair distance: {0} (should be around 0.276 nm)".format(dmin))
                self.density = density0 / (0.276 / dmin)**3
                # self.density = density0
        else:
            self.density = density

        self.logger.info("Target Density: {0}".format(self.density))
        self.logger.info("Original Density: {0}".format(density0))

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density)**(1.0 / 3.0)
        self.cell.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        self.logger.info("Bond length (scaled, nm): {0}".format(self.bondlen))

        # ================================================================
        # double_network: True or False
        #   This is a special option for ices VI and VII that have
        #   interpenetrating double network.
        #   GenIce's fast depolarization algorithm fails in some case.
        #
        try:
            self.double_network = lat.double_network
        except AttributeError:
            self.double_network = False

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
            self.logger.info("Orientations of some edges are fixed.")
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

        
    def analyze_ice(self, water_type, formatter, noise=0.):
        """
        Protocol for analice
        """
        hooks = formatter.hooks
        arg   = formatter.arg

        maxstage = max(0, *hooks.keys())

        if 0 in hooks:
            hooks[0](self, arg)
        elif arg != "":
            logger.info("Arguments are given but the module does not accept them.")
            if "usage" in formatter.__dict__:
                formatter.usage()
            else:
                for line in formatter.__doc__.splitlines():
                    logger.info("  "+line)

        self.stage1(noise)

        if 1 in hooks:
            abort = hooks[1](self)
            if maxstage < 2 or abort:
                return

        if self.rotmatrices is None:
            res = self.stage2()

        if 2 in hooks:
            abort = hooks[2](self)
            if maxstage < 3 or abort:
                return

        if self.rotmatrices is None:
            self.stage3()

        if 3 in hooks:
            abort = hooks[3](self)
            if maxstage < 4 or abort:
                return

        self.stage4()

        if 4 in hooks:
            abort = hooks[4](self)
            if maxstage < 5 or abort:
                return

        # molecular orientation should be given in the loader.
        if self.rotmatrices is None:
            self.stage5()

        if 5 in hooks:
            abort = hooks[5](self)
            if maxstage < 6 or abort:
                return

        self.stage6(water_type)

        if 6 in hooks:
            abort = hooks[6](self)
            if maxstage < 7 or abort:
                return

        # self.stage7_analice(guests)
        if 7 in hooks:
            hooks[7](self)

    def stage1(self,
               noise=0.):
        """
        Do nothing.

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        """

        self.logger.info("Stage1: (...)")
        self.reppositions = self.waters

        # This must be done before the replication of the cell.
        self.logger.info("  Number of water molecules: {0}".format(len(self.reppositions)))

        # self.graph = self.prepare_random_graph(self.fixed)
        self.graph = self.prepare_random_graph(self.pairs)

        # scale the cell
        self.repcell = Cell(self.cell.mat)

        # self.repcell.scale2(self.rep)
        # add small perturbations to the molecular positions.
        if noise > 0.0:
            self.logger.info("  Add noise: {0}.".format(noise))
            perturb = np.random.normal(loc=0.0,
                                       scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                                       size=self.reppositions.shape)
            self.reppositions += self.repcell.abs2rel(perturb)

        self.logger.info("Stage1: end.")

    def stage4(self):
        """
        Depolarize.

        Provided variables:
        spacegraph: depolarized network with node positions.
        yapresult:  Animation of the depolarization process in YaPlot format.
        """

        self.logger.info("Stage4: (...)")
        self.yapresult = ""
        self.spacegraph = dg.SpaceIceGraph(self.graph,
                                           coord=self.reppositions,
                                           ignores=self.graph.ignores)
        self.logger.info("Stage4: end.")
