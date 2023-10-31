"""
GenIce class
"""

import itertools as it
from collections import defaultdict
from logging import getLogger
from typing import Type, Union

import numpy as np
import pairlist as pl
import networkx as nx

from genice2 import cage

# from genice2 import digraph_unused as dg
from genice2.cell import Cell, rel_wrap
from genice2.decorators import banner, timeit
from genice2.formats import Format
from genice2.lattices import Lattice
import genice_core

# A virtual monatomic molecule
from genice2.molecules import Molecule, arrange, monatom, one
from genice2.plugin import Group, safe_import
from genice2.valueparser import (
    parse_cages,
    parse_pairs,
    plugin_option_parser,
    put_in_array,
)


def assume_tetrahedral_vectors(v):
    """
    Assume missing vectors at a tetrahedral node.

    Given: known vectors.
    Returns: assumed vectors
    """

    assert len(v) > 0

    if len(v) == 3:
        return [-(v[0] + v[1] + v[2])]

    if len(v) == 2:
        y = v[1] - v[0]
        y /= np.linalg.norm(y)
        z = v[1] + v[0]
        z /= np.linalg.norm(z)
        x = np.cross(y, z)
        v2 = (x * 8.0**0.5 - z) / 3.0
        v3 = (-x * 8.0**0.5 - z) / 3.0
        return [v2, v3]

    if len(v) == 1:
        vr = np.array([random.random() for i in range(3)])
        vr /= np.linalg.norm(vr)
        z = v[0] / np.linalg.norm(v[0])
        x = np.cross(z, vr)
        y = np.cross(z, x)
        x1 = -x / 2 + 3.0**0.5 * y / 2
        x2 = -x / 2 - 3.0**0.5 * y / 2
        return [x, x1, x2]

    return []


def ice_rule(dg: nx.DiGraph, strict=False) -> bool:
    if strict:
        for node in dg:
            if dg.in_degree(node) != 2 or dg.out_degree(node) != 2:
                return False
    else:
        for node in dg:
            if dg.in_degree(node) > 2 or dg.out_degree(node) > 2:
                return False
    return True


def orientations(coord, graph, cell, immutables: set):
    """
    Does not work when two OHs are colinear
    """

    logger = getLogger()
    # just for a test of pure water
    assert len(coord) == graph.number_of_nodes(), (len(coord), graph.number_of_nodes())

    # 通常の氷であればアルゴリズムを高速化できる。

    if ice_rule(graph, strict=True) and len(immutables) == 0:
        # fast track
        rotmatrices = np.zeros([len(list(graph)), 3, 3])

        neis = np.zeros([len(list(graph)), 2], dtype=int)
        for node in graph:
            neis[node] = list(graph.successors(node))
        # array of donating vectors
        v0 = coord[neis[:, 0]] - coord[:]
        v0 -= np.floor(v0 + 0.5)
        v0 = v0 @ cell.mat
        v0 /= np.linalg.norm(v0, axis=1)[:, np.newaxis]
        v1 = coord[neis[:, 1]] - coord[:]
        v1 -= np.floor(v1 + 0.5)
        v1 = v1 @ cell.mat
        v1 /= np.linalg.norm(v1, axis=1)[:, np.newaxis]
        # intramolecular axes
        y = v1 - v0
        y /= np.linalg.norm(y, axis=1)[:, np.newaxis]
        z = v0 + v1
        z /= np.linalg.norm(z, axis=1)[:, np.newaxis]
        x = np.cross(y, z, axisa=1, axisb=1)
        rotmatrices[:, 0, :] = x
        rotmatrices[:, 1, :] = y
        rotmatrices[:, 2, :] = z
        return rotmatrices

    else:
        rotmatrices = []
        for node in range(graph.number_of_nodes()):
            if node in immutables:
                # for dopants; do not rotate
                rotmat = np.identity(3)
            else:
                vsucc = [
                    cell.rel2abs(rel_wrap(coord[x] - coord[node]))
                    for x in graph.successors(node)
                ]

                if len(vsucc) < 2:  # TSL
                    vpred = [
                        cell.rel2abs(rel_wrap(coord[x] - coord[node]))
                        for x in graph.predecessors(node)
                    ]
                    vsucc = [x / np.linalg.norm(x) for x in vsucc]
                    vpred = [x / np.linalg.norm(x) for x in vpred]

                    if len(vpred) > 2:
                        # number of incoming bonds should be <= 2
                        vpred = vpred[:2]
                    vcomp = assume_tetrahedral_vectors(vpred + vsucc)
                    logger.debug(
                        f"Node {node} vcomp {vcomp} vsucc {vsucc} vpred {vpred}"
                    )
                    vsucc = (vsucc + vcomp)[:2]

                logger.debug(f"Node {node} vsucc {vsucc}")
                assert 2 <= len(vsucc), "Probably a wrong ice network."
                # normalize vsucc
                vsucc[0] /= np.linalg.norm(vsucc[0])
                vsucc[1] /= np.linalg.norm(vsucc[1])
                y = vsucc[1] - vsucc[0]
                y /= np.linalg.norm(y)
                z = (vsucc[0] + vsucc[1]) / 2
                z /= np.linalg.norm(z)
                x = np.cross(y, z)
                # orthogonality check
                # logger.debug((x@x,y@y,z@z,x@y,y@z,z@x))
                rotmat = np.vstack([x, y, z])

            rotmatrices.append(rotmat)

    return rotmatrices


def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99

    if pairs is None:
        iter = it.combinations(coord, 2)
    else:
        iter = [(coord[i], coord[j]) for i, j in pairs]

    for c1, c2 in iter:
        r = cell.rel2abs(rel_wrap(c1 - c2))
        rr = r @ r

        if rr < dmin:
            dmin = rr

    return dmin**0.5


def replicate_groups(groups, waters, cagepos, rep):
    """
    This is not that easy.
    """
    logger = getLogger()
    # Storage for replicated groups
    newgroups = defaultdict(dict)

    for root, cages in groups.items():
        # Position of root (water) (fractional)
        root_pos = waters[root]

        for cage, group_name in cages.items():
            # Position of the cage (fractional)
            cage_pos = cagepos[cage]
            # Relative position of the cage
            delta = rel_wrap(cage_pos - root_pos)
            # (Image) cell that the cage resides
            gcell = np.floor(root_pos + delta)

            for x in range(rep[0]):
                for y in range(rep[1]):
                    for z in range(rep[2]):
                        r = np.array((x, y, z))
                        # label of the root (water) in the replica
                        newroot = root + len(waters) * (x + rep[0] * (y + rep[1] * z))
                        # replicated cell in which the cage resides.
                        # modulo by positive number is always positive.
                        cr = (r + gcell) % rep
                        newcage = cage + len(cagepos) * (
                            cr[0] + rep[0] * (cr[1] + rep[1] * cr[2])
                        )
                        newcage = int(newcage)
                        newgroups[newroot][newcage] = group_name
                        # logger.info(("root",newroot,"newcage", newcage))
    return newgroups


def replicate_labeldict(labels, nmol, rep):
    newlabels = dict()

    for j in labels:
        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    newj = j + nmol * (x + rep[0] * (y + rep[1] * z))
                    newlabels[newj] = labels[j]

    return newlabels


def replicate_labels(labels, nmol, rep):
    newlabels = set()

    for j in labels:
        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    newj = j + nmol * (x + rep[0] * (y + rep[1] * z))
                    newlabels.add(newj)

    return newlabels


@timeit
def replicate_graph(graph, positions, rep, fixed: list):
    # repgraph = dg.IceGraph()
    repgraph = nx.Graph()
    fixed = nx.DiGraph(fixed)
    fixedEdges = nx.DiGraph()
    nmol = positions.shape[0]

    for i, j in graph.edges(data=False):
        if fixed.has_edge(i, j):
            edge_fixed = True
        elif fixed.has_edge(j, i):
            edge_fixed = True
            i, j = j, i
        else:
            edge_fixed = False
        # positions in the unreplicated cell
        vec = positions[j] - positions[i]
        delta = np.floor(vec + 0.5).astype(int)

        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    xi = (x + delta[0]) % rep[0]
                    yi = (y + delta[1]) % rep[1]
                    zi = (z + delta[2]) % rep[2]
                    newi = i + nmol * (xi + rep[0] * (yi + rep[1] * zi))
                    newj = j + nmol * (x + rep[0] * (y + rep[1] * z))

                    if edge_fixed:  # fixed pair
                        fixedEdges.add_edge(newi, newj)
                    repgraph.add_edge(newi, newj)

    return repgraph, fixedEdges


def replicate_positions(positions, rep):
    repx = positions.copy()

    for m in range(1, rep[0]):
        v = np.array([m, 0, 0])
        repx = np.concatenate((repx, positions + v))

    repy = repx.copy()

    for m in range(1, rep[1]):
        v = np.array([0, m, 0])
        repy = np.concatenate((repy, repx + v))

    repz = repy.copy()

    for m in range(1, rep[2]):
        v = np.array([0, 0, m])
        repz = np.concatenate((repz, repy + v))

    return repz / rep


def neighbor_cages_of_dopants(dopants, waters, cagepos, cell):
    """
    Just shows the environments of the dopants
    """
    # logger = getLogger()
    dnei = defaultdict(set)

    for site, name in dopants.items():
        org = waters[site]

        for i, pos in enumerate(cagepos):
            # Displacement (relative)
            a = cell.rel2abs(rel_wrap(pos - org))
            sqdistance = a @ a

            if sqdistance < 0.57**2:
                dnei[site].add(i)
                # logger.info((i,cagepos[i]))

    return dnei


class GenIce:
    """
    The core of GenIce.

    lat:     An instance of a Lattice class.
    density: Density of target ice in g/cm3.
    rep:     Repetition of the cell. A tuple of three integers.
    anions, cations:
             The locations of monovalent anions and cations that replace
             the water molecules.
    spot_guests:
             The locations of guest molecules that occupy the specified cages.
    spot_group: (EXPERIMENTAL)
             The locations of functional groups that occupy the cages.
    as_is:   Avoids shuffling of the orientations of water molecules.
    signature: A text that is inserted in the output.
    """

    @timeit
    @banner
    def __init__(
        self,
        lat: Type[Lattice],
        signature: str = "",
        density: float = 0,
        rep=(1, 1, 1),
        cations: dict = {},
        anions: dict = {},
        spot_guests: dict = {},
        spot_groups: dict = {},
        asis: bool = False,
        shift=(0.0, 0.0, 0.0),
        # seed=1000,
    ):
        """
        Constructor of GenIce.

        Arguments:
            lat:        The ice lattice.
            signature:  A string for a signature.
            density:    Target density.
            rep:        Cell repetitions.
            cations:    Labels of water molecules that are replaced by the cations.
            anions:
            spot_guests:Labels of cages in which a guest is placed.
            spot_groups:Labels of cages in which a group is placed.
            asis:       Do not modify the orientations of the hydrogen bonds.
            shift:      A fractional value to be added to the positions.

        """

        # このconstructorが大きすぎ、複雑すぎ。
        # self変数も多すぎ。
        logger = getLogger()
        self.rep = rep
        self.asis = asis
        self.cations = cations
        self.anions = anions
        self.spot_guests = spot_guests
        self.spot_groups = spot_groups

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

        # 変数名に1がついているのはunit cell

        # ================================================================
        # cell: cell dimension
        #   see parse_cell for syntax.
        #
        self.cell1 = Cell(lat.cell)

        # ================================================================
        # waters: positions of water molecules
        #
        self.waters1 = put_in_array(lat.waters)
        logger.debug(f"Waters: {len(self.waters1)}")
        self.waters1 = self.waters1.reshape((-1, 3))

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative", i.e. in fractional coordinates
        #
        if lat.coord == "absolute":
            self.waters1 = self.cell1.abs2rel(self.waters1)

        # shift of the origin
        self.waters1 = np.array(self.waters1) + np.array(shift)
        # fractional coordinate between [0, 1)
        self.waters1 -= np.floor(self.waters1)

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
            logger.info(f"Bond length (specified): {self.bondlen}")
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            # very rough estimate of the pair distances assuming that the particles distribute homogeneously.
            rc = (volume / nmol) ** (1 / 3) * 1.5
            p = pl.pairs_iter(
                self.waters1, maxdist=rc, cell=self.cell1.mat, distance=False
            )
            self.bondlen = 1.1 * shortest_distance(self.waters1, self.cell1, pairs=p)
            logger.info(f"Bond length (estim.): {self.bondlen}")

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
                    f"Closest pair distance: {dmin} (should be around 0.276 nm)"
                )
                self.density = density0 / (0.276 / dmin) ** 3
                # self.density = density0
        else:
            self.density = density

        logger.info(f"Target Density: {self.density}")
        logger.info(f"Original Density: {density0}")

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density) ** (1.0 / 3.0)
        self.cell1.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        logger.info(f"Bond length (scaled, nm): {self.bondlen}")

        # ================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        self.cagepos1 = None
        self.cagetype1 = None

        if "cages" in lat.__dict__:
            self.cagepos1, self.cagetype1 = parse_cages(lat.cages)
            logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            # pre-parsed data
            self.cagepos1, self.cagetype1 = np.array(lat.cagepos), lat.cagetype
        if self.cagepos1 is not None:
            self.cagepos1 = np.array(self.cagepos1) + np.array(shift)
            self.cagepos1 -= np.floor(self.cagepos1)

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

        if "dopeIonsToUnitCell" in lat.__dict__:
            self.dopeIonsToUnitCell = lat.dopeIonsToUnitCell
        else:
            self.dopeIonsToUnitCell = None
        self.dopants1 = set()

        self.immutables1 = set(self.dopants1)

        # if asis, make pairs to be fixed.
        if self.asis and len(self.fixed1) == 0:
            self.fixed1 = self.pairs1

        # filled cages
        self.filled_cages = set()

        # groups info
        self.groups1 = defaultdict(dict)

    def generate_ice(
        self,
        formatter: Type[Format],
        water: Union[Type[Molecule], None] = None,
        guests={},
        depol="strict",
        noise=0.0,
        assess_cages=False,
    ):
        """
        Generate an ice structure and dump it with the aid of a formatter plugin.

            formatter: genice2.format.Format() class
            water:     genice2.molecules.Molecule() class
            assess_cages:   Cages will be assessed on the fly instead of
                        pre-specified in the lattice plugin.
        """

        logger = getLogger()

        # in old syntax, the arguments water and formatter were mandatory, but
        # in new syntax, water is optional and their order is exchanged.
        # therefore i prepare a backward compatibility.
        from genice2.molecules import Molecule

        if isinstance(formatter, Molecule):
            formatter, water = water, formatter
            logger.warn(
                "generate_ice(water, formatter) is deprecated. "
                "New syntax is: generate_ice(formatter, water=water)."
            )

        def Stages():
            hooks = formatter.hooks()
            maxstage = max(0, *hooks.keys())

            if 0 in hooks:
                abort = hooks[0](self)
                if maxstage < 1 or abort:
                    return

            self.Stage1(noise, assess_cages=assess_cages)

            if 1 in hooks:
                abort = hooks[1](self)
                if maxstage < 2 or abort:
                    return

            self.Stage2()

            # Count bonds
            nfixed = self.fixedEdges.number_of_edges()
            num_hb_disorder = self.graph.number_of_edges() - nfixed
            logger.info(f"  Number of pre-oriented hydrogen bonds: {nfixed}")
            logger.info(f"  Number of unoriented hydrogen bonds: {num_hb_disorder}")
            logger.info(
                "  Number of hydrogen bonds: {0} (regular num: {1})".format(
                    nfixed + num_hb_disorder, len(self.reppositions) * 2
                )
            )

            if 2 in hooks:
                abort = hooks[2](self)
                if maxstage < 3 or abort:
                    return

            # GenIce-core
            self.Stage34E(depol=depol)

            if 3 in hooks:
                abort = hooks[3](self)
                if maxstage < 4 or abort:
                    return

            if 4 in hooks:
                abort = hooks[4](self)
                if maxstage < 5 or abort:
                    return

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

            self.Stage7(guests)

            if 7 in hooks:
                hooks[7](self)

        abort = Stages()
        if not abort:
            return formatter.dump()

    @timeit
    @banner
    def Stage1(self, noise=0.0, assess_cages=False):
        """
        Replicate water molecules to make a repeated cell.

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        repcagetype: replicated cage types array
        repcagepos:  replicated cage positions (CoM, relative)
        cagetypes:   set of cage types
        """

        logger = getLogger()
        self.reppositions = replicate_positions(self.waters1, self.rep)

        # This must be done before the replication of the cell.
        logger.info("  Number of water molecules: {0}".format(len(self.reppositions)))

        if self.pairs1 is None:
            self.pairs1 = self.prepare_pairs()

        self.graph1 = nx.Graph(self.pairs1)
        logger.info(f"  Number of water nodes: {self.graph1.number_of_nodes()}")

        # scale the cell
        self.repcell = Cell(self.cell1.mat)
        self.repcell.scale2(self.rep)

        if noise > 0.0:
            logger.info(f"  Add noise: {noise}.")
            perturb = np.random.normal(
                loc=0.0,
                scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                size=self.reppositions.shape,
            )
            self.reppositions += self.repcell.abs2rel(perturb)

        if assess_cages:
            logger.info("  Assessing the cages...")
            self.cagepos1, self.cagetype1 = cage.assess_cages(self.graph1, self.waters1)
            logger.info("  Done assessment.")

        if self.cagepos1 is not None and self.cagepos1.shape[0] > 0:
            logger.info("  Hints:")
            self.repcagepos = replicate_positions(self.cagepos1, self.rep)
            nrepcages = self.repcagepos.shape[0]
            self.repcagetype = [
                self.cagetype1[i % len(self.cagetype1)] for i in range(nrepcages)
            ]
            self.cagetypes = defaultdict(set)

            for i, typ in enumerate(self.repcagetype):
                self.cagetypes[typ].add(i)

            # INFO for cage types
            logger.info("    Cage types: {0}".format(list(self.cagetypes)))

            for typ, cages in self.cagetypes.items():
                logger.info(f"    Cage type {typ}: {cages}")
            # Up here move to stage 1.
        else:
            self.repcagetype = None
            self.repcagepos = None
            self.cagetypes = None

    @timeit
    @banner
    def Stage2(self):
        """
        Make a random graph and replicate.

        Provided variables:
        dopants:
        groups:  replicated positions of the chemical groups (CoM, relative)
        filled_cages:
        graph:   replicated network topology (bond orientation may be random)
        """

        logger = getLogger()

        # Some edges are directed when ions are doped.
        if self.dopeIonsToUnitCell is not None:
            self.dopeIonsToUnitCell(self)  # may be defined in the plugin

        # Replicate the dopants in the unit cell
        self.dopants = replicate_labeldict(self.dopants1, len(self.waters1), self.rep)
        self.groups = replicate_groups(
            self.groups1, self.waters1, self.cagepos1, self.rep
        )

        # self.groups_info(self.groups)
        for _, cages in self.groups.items():
            self.filled_cages |= set(cages)

        # logger.info(("filled",self.filled_cages))
        # Replicate the graph
        self.graph, self.fixedEdges = replicate_graph(
            self.graph1, self.waters1, self.rep, self.fixed1
        )
        # replicate "immutables" == dopants list in the graph
        self.repimmutables = replicate_labels(
            self.immutables1, self.waters1.shape[0], self.rep
        )

        # Dope ions by options.
        if len(self.anions) > 0:
            logger.info(f"  Anionize: {self.anions}.")

            for site, name in self.anions.items():
                # self.graph.anionize(site)
                for nei in self.graph[site]:
                    self.fixedEdges.add_edge(nei, site)
                self.dopants[site] = name

        if len(self.cations) > 0:
            logger.info(f"  Cationize: {self.cations}.")

            for site, name in self.cations.items():
                # self.graph.cationize(site)
                for nei in self.graph[site]:
                    self.fixedEdges.add_edge(site, nei)
                self.dopants[site] = name

    @timeit
    @banner
    def Stage34E(self, depol: str = "none"):
        """
        Make a directed ice graph from an undirected ice graph.
        """

        logger = getLogger()

        if depol == "none":
            iter = 0
        else:
            iter = 1000
        d = genice_core.ice_graph(
            self.graph.to_undirected(),
            vertexPositions=self.reppositions,
            isPeriodicBoundary=True,
            dipoleOptimizationCycles=iter,
            fixedEdges=self.fixedEdges,
        )

        self.digraph = nx.compose(d, self.fixedEdges)

        logger.debug("  Depolarized in Stage34E.")

    @timeit
    @banner
    def Stage5(self):
        """
        Prepare orientations for rigid molecules.

        Provided variables:
        reppositions: molecular positions.
        rotmatrices:  rotation matrices for water molecules
        """

        logger = getLogger()

        # determine the orientations of the water molecules based on edge
        # directions.
        self.rotmatrices = orientations(
            self.reppositions, self.digraph, self.repcell, self.repimmutables
        )

        # Activate it.
        # logger.info("The network is not specified.  Water molecules will be orinented randomly.")
        # rotmatrices = [rigid.rand_rotation_matrix() for pos in positions]

    @timeit
    @banner
    def Stage6(self, water):
        """
        Arrange atoms of water and replacements

        Provided variables:
        atoms: atomic positions of water molecules. (absolute)
        """

        logger = getLogger()

        # assert audit_name(water_type), "Dubious water name: {0}".format(water_type)
        # water = importlib.import_module("genice.molecules."+water_type)
        # water = safe_import("molecule", water_type)

        try:
            mdoc = water.__doc__.splitlines()
        except BaseException:
            mdoc = []

        for line in mdoc:
            logger.info("  " + line)

        self.universe = []
        self.universe.append(
            arrange(
                self.reppositions,
                self.repcell,
                self.rotmatrices,
                water,
                immutables=set(self.dopants),
            )
        )

    @timeit
    @banner
    def Stage7(self, guests):
        """
        Arrange guest atoms.

        Provided variables:
        atoms: atomic positions of all molecules.
        """

        logger = getLogger()

        if self.cagepos1 is not None:
            # the cages around the dopants.
            dopants_neighbors = dopants_info(
                self.dopants, self.reppositions, self.repcagepos, self.repcell
            )

            # put the (one-off) groups
            if len(self.spot_groups) > 0:
                # process the -H option
                for cage, group_to in self.spot_groups.items():
                    group, root = group_to.split(":")
                    self.add_group(cage, group, int(root))

            molecules = defaultdict(list)

            if len(self.spot_guests) > 0:
                # process the -G option
                for cage, molec in self.spot_guests.items():
                    molecules[molec].append(cage)
                    self.filled_cages.add(cage)

            # process the -g option
            for cagetype, contents in guests.items():
                assert cagetype in self.cagetypes, f"Nonexistent cage type: {cagetype}"
                resident = dict()
                rooms = list(self.cagetypes[cagetype] - self.filled_cages)

                for room in rooms:
                    resident[room] = None

                vacant = len(rooms)

                for molec, frac in contents.items():
                    nmolec = int(frac * len(rooms) + 0.5)
                    vacant -= nmolec
                    assert vacant >= 0, "Too many guests."
                    remain = nmolec
                    movedin = []

                    while remain > 0:
                        r = np.random.randint(0, len(rooms))
                        room = rooms[r]

                        if resident[room] is None:
                            resident[room] = molec
                            molecules[molec].append(room)
                            movedin.append(room)
                            remain -= 1

            # Now ge got the address book of the molecules.
            if len(molecules):
                logger.info("  Summary of guest placements:")
                guests_info(self.cagetypes, molecules)

            if len(self.spot_groups) > 0:
                logger.info("  Summary of groups:")
                groups_info(self.groups)

            # semi-guests
            for root, cages in self.groups.items():
                assert root in self.dopants
                name = self.dopants[root]
                molname = f"G{root}"
                pos = self.reppositions[root]
                rot = self.rotmatrices[root]
                self.universe.append(monatom(pos, self.repcell, name))
                del self.dopants[root]  # processed.
                logger.debug((root, cages, name, molname, pos, rot))

                for cage, group in cages.items():
                    # assert group in self.groups_placer
                    assert cage in dopants_neighbors[root]
                    cage_center = self.repcagepos[cage]
                    self.universe.append(
                        Group(group).arrange_atoms(
                            cage_center, pos, self.repcell, molname, origin_atom=name
                        )
                    )

            # molecular guests
            for molec, cages in molecules.items():
                guest_type, guest_options = plugin_option_parser(molec)
                logger.debug(f"Guest type: {guest_type}")
                gmol = safe_import("molecule", guest_type).Molecule(**guest_options)

                try:
                    mdoc = gmol.__doc__.splitlines()
                except BaseException:
                    mdoc = []
                for line in mdoc:
                    logger.info("  " + line)
                cage_center = [self.repcagepos[i] for i in cages]
                cmat = [np.identity(3) for i in cages]
                self.universe.append(arrange(cage_center, self.repcell, cmat, gmol))

        # Assume the dopant is monatomic and replaces one water molecule
        atomset = defaultdict(set)
        for label, name in self.dopants.items():
            atomset[name].add(label)

        for name, labels in atomset.items():
            pos = [self.reppositions[i] for i in sorted(labels)]
            rot = [self.rotmatrices[i] for i in sorted(labels)]
            oneatom = one.Molecule(label=name)
            self.universe.append(arrange(pos, self.repcell, rot, oneatom))

    def prepare_pairs(self):
        logger = getLogger()

        logger.info("  Pairs are not given explicitly.")
        logger.info("  Estimating the bonds according to the pair distances.")

        logger.debug(f"Bondlen: {self.bondlen}")
        # make bonded pairs according to the pair distance.
        # make before replicating them.
        return [
            (i, j)
            for i, j in pl.pairs_iter(
                self.waters1, self.bondlen, self.cell1.mat, distance=False
            )
        ]

    def add_group(self, cage, group, root):
        self.groups[root][cage] = group
        self.filled_cages.add(cage)

    def __del__(self):
        logger = getLogger()
        logger.info("Completed.")


def dopants_info(dopants=None, waters=None, cagepos=None, cell=None):
    logger = getLogger()
    if dopants is None:
        dopants = self.dopants

    if waters is None:
        waters = self.waters

    if cagepos is None:
        cagepos = self.cagepos

    if cell is None:
        cell = self.cell

    dopants_neighbors = neighbor_cages_of_dopants(dopants, waters, cagepos, cell)

    for dopant, cages in dopants_neighbors.items():
        logger.info(f"    Cages adjacent to dopant {dopant}: {cages}")

    return dopants_neighbors


def groups_info(groups):
    logger = getLogger()
    for root, cages in groups.items():
        for cage, group in cages.items():
            logger.info(f"    Group {group} of dopant {root} in cage {cage}")


def guests_info(cagetypes, molecules):
    logger = getLogger()
    for cagetype, cageid in cagetypes.items():
        logger.info(f"    Guests in cage type {cagetype}:")

        for molec, cages in molecules.items():
            cages = set(cages)
            cages &= cageid

            if len(cages):
                logger.info(f"      {molec} * {len(cages)} @ {cages}")
