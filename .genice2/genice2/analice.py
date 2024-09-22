from collections import defaultdict
from logging import getLogger

import numpy as np
import pairlist as pl
import networkx as nx

# from genice2 import digraph_unused as dg
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

        # unit cellに関する変数には1をつける。

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
        self.waters1 = put_in_array(lat.waters)
        logger.debug("Waters: {0}".format(len(self.waters1)))
        self.waters1 = self.waters1.reshape((self.waters1.size // 3, 3))

        # ================================================================
        # cell: cell dimension
        #   see parse_cell for syntax.
        #
        self.cell1 = Cell(lat.cell)

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative"
        #
        if lat.coord == "absolute":
            self.waters1 = self.cell1.abs2rel(self.waters1)

        self.waters1 = np.array([w - np.floor(w) for w in self.waters1])

        # ================================================================
        # pairs: specify the pairs of molecules that are connected.
        #   Bond orientation will be shuffled later
        #   unless it is "fixed".
        #
        self.pairs1 = None

        try:
            self.pairs1 = parse_pairs(lat.pairs)
        except AttributeError:
            logger.info("HB connectivity is not defined.")

        # ================================================================
        # bondlen: specify the bond length threshold.
        #   This is used when "pairs" are not specified.
        #   It is applied to the original positions of molecules (before density setting).
        #
        nmol = self.waters1.shape[0]  # nmol in a unit cell
        volume = self.cell1.volume()  # volume of a unit cell in nm**3
        self.bondlen = None

        try:
            self.bondlen = lat.bondlen
            logger.info("Bond length (specified): {0}".format(self.bondlen))
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            # assume that the particles distribute homogeneously.
            rc = (volume / nmol) ** (1 / 3) * 1.5
            p = pl.pairs_iter(
                self.waters1, maxdist=rc, cell=self.cell1.mat, distance=False
            )
            self.bondlen = 1.1 * shortest_distance(self.waters1, self.cell1, pairs=p)
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
                    "Density is not specified. Assume the density from lattice."
                )
                dmin = shortest_distance(self.waters1, self.cell1)
                logger.info(
                    "Closest pair distance: {0} (should be around 0.276 nm)".format(
                        dmin
                    )
                )
                self.density = density0 / (0.276 / dmin) ** 3
                # self.density = density0
        else:
            self.density = density

        logger.info("Target Density: {0}".format(self.density))
        logger.info("Original Density: {0}".format(density0))

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density) ** (1.0 / 3.0)
        self.cell1.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        logger.info("Bond length (scaled, nm): {0}".format(self.bondlen))

        # ================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        self.cagepos1 = None
        self.cagetype1 = None

        # ================================================================
        # fixed: specify the bonds whose directions are fixed.
        #   you can specify them in pairs at a time.
        #   You can also leave it undefined.
        #
        self.fixed1 = []
        try:
            self.fixed1 = parse_pairs(lat.fixed)
            logger.info("Orientations of some edges are fixed.")
        except AttributeError:
            pass

        self.dopeIonsToUnitCell = None
        self.dopants1 = set()
        # analiceではセルの複製はしない?
        self.dopants = self.dopants1

        # if asis, make pairs to be fixed.
        if self.asis and len(self.fixed1) == 0:
            self.fixed1 = self.pairs1

        # analiceではセルの複製はしない?
        self.fixedEdges = nx.DiGraph(self.fixed1)

        # filled cages
        self.filled_cages1 = set()

        # groups info
        self.groups1 = defaultdict(dict)

        # unused but required
        self.anions1 = []
        self.cations1 = []

    def analyze_ice(self, water, formatter, noise=0.0):
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

            # if self.rotmatrices is None:
            #     self.Stage3()

            # if 3 in hooks:
            #     abort = hooks[3](self)
            #     if maxstage < 4 or abort:
            #         return

            # self.Stage4()

            # if 4 in hooks:
            #     abort = hooks[4](self)
            #     if maxstage < 5 or abort:
            #         return

            # GenIce-core
            self.Stage34E()

            if 3 in hooks:
                abort = hooks[3](self)
                if maxstage < 4 or abort:
                    return

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
    def Stage1(self, noise=0.0):
        """Preparation.

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        """

        logger = getLogger()
        self.reppositions = self.waters1

        # This must be done before the replication of the cell.
        logger.info("  Number of water molecules: {0}".format(len(self.reppositions)))

        # self.graph = self.prepare_random_graph(self.fixed)
        if self.pairs1 is None:
            self.pairs1 = self.prepare_pairs()

        self.graph1 = nx.Graph(self.pairs1)
        logger.info(f"  Number of water nodes: {self.graph1.number_of_nodes()}")
        self.graph = self.graph1

        # scale the cell
        self.repcell = Cell(self.cell1.mat)

        # self.repcell.scale2(self.rep)
        # add small perturbations to the molecular positions.
        if noise > 0.0:
            logger.info("  Add noise: {0}.".format(noise))
            perturb = np.random.normal(
                loc=0.0,
                scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                size=self.reppositions.shape,
            )
            self.reppositions += self.repcell.abs2rel(perturb)

    @timeit
    @banner
    def Stage4(self):
        """Make a spacegraph.

        Provided variables:
        spacegraph: depolarized network with node positions.
        yapresult:  Animation of the depolarization process in YaPlot format.
        """

        self.yapresult = ""
        self.digraph = self.graph
