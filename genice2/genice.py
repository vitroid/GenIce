#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GenIce class
"""

import itertools as it
import random
from collections import defaultdict
from logging import getLogger
from typing import Type, Union

import numpy as np
import pairlist as pl
import tilecycles as tc

from genice2 import cage
from genice2 import digraph as dg
from genice2.cell import Cell, rel_wrap
from genice2.decorators import banner, timeit
from genice2.formats import Format
from genice2.lattices import Lattice
from genice2.molecules import Molecule, arrange, monatom, one
from genice2.plugin import Group, safe_import
from genice2.valueparser import (flatten, parse_cages, parse_pairs,
                                 plugin_option_parser, put_in_array)


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


def orientations(coord, graph, cell):
    """
    Does not work when two OHs are colinear
    """

    logger = getLogger()
    # just for a test of pure water
    assert len(coord) == graph.number_of_nodes()

    # 通常の氷であればアルゴリズムを高速化できる。

    if graph.isZ22() and len(graph.immutables) == 0:
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
            if node in graph.immutables:
                # for dopants; do not rotate
                rotmat = np.identity(3)
            else:
                vsucc = [cell.rel2abs(rel_wrap(coord[x] - coord[node]))
                         for x in graph.successors(node)]

                if len(vsucc) < 2:  # TSL
                    vpred = [cell.rel2abs(rel_wrap(coord[x] - coord[node]))
                             for x in graph.predecessors(node)]
                    vsucc = [x / np.linalg.norm(x) for x in vsucc]
                    vpred = [x / np.linalg.norm(x) for x in vpred]

                    if len(vpred) > 2:
                        # number of incoming bonds should be <= 2
                        vpred = vpred[:2]
                    vcomp = assume_tetrahedral_vectors(vpred + vsucc)
                    logger.debug(
                        f"Node {node} vcomp {vcomp} vsucc {vsucc} vpred {vpred}")
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
                        newroot = root + len(waters) * \
                            (x + rep[0] * (y + rep[1] * z))
                        # replicated cell in which the cage resides.
                        # modulo by positive number is always positive.
                        cr = (r + gcell) % rep
                        newcage = cage + \
                            len(cagepos) * (cr[0] + rep[0]
                                            * (cr[1] + rep[1] * cr[2]))
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
def replicate_graph(graph, positions, rep):
    repgraph = dg.IceGraph()
    nmol = positions.shape[0]

    for i, j in graph.edges(data=False):
        # positions in the unreplicated cell
        vec = positions[j] - positions[i]
        delta = np.floor(vec + 0.5).astype(int)
        edge_fixed = graph[i][j]['fixed']

        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    xi = (x + delta[0]) % rep[0]
                    yi = (y + delta[1]) % rep[1]
                    zi = (z + delta[2]) % rep[2]
                    newi = i + nmol * (xi + rep[0] * (yi + rep[1] * zi))
                    newj = j + nmol * (x + rep[0] * (y + rep[1] * z))

                    if edge_fixed:  # fixed pair
                        repgraph.add_edge(newi, newj, fixed=True)
                    else:

                        # shuffle the bond directions
                        if 0 == random.randint(0, 1):
                            repgraph.add_edge(newi, newj, fixed=False)
                        else:
                            repgraph.add_edge(newj, newi, fixed=False)

    # replicate "immutables" == dopants list in the graph
    repgraph.immutables = replicate_labels(graph.immutables, nmol, rep)

    return repgraph


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
    #logger = getLogger()
    dnei = defaultdict(set)

    for site, name in dopants.items():
        org = waters[site]

        for i, pos in enumerate(cagepos):
            #Displacement (relative)
            a = cell.rel2abs(rel_wrap(pos - org))
            sqdistance = a @ a

            if sqdistance < 0.57**2:
                dnei[site].add(i)
                # logger.info((i,cagepos[i]))

    return dnei


# They should be separate plugins in the future.


# def pentyl(cage_center, root_position, cell, molname):
#     """
#     put a butyl group rooted at root_position toward cage_center.
#     """
#     return Alkyl(
#         cage_center, root_position, cell, molname, [
#             "Ma", [
#                 "Mb", [
#                     "Mc", [
#                         "Md", "Me"]]]])


# def propyl(cage_center, root_position, cell, molname):
#     return Alkyl(cage_center, root_position, cell, molname, ["Ma", ["Mb", "Mc"]])


# def ethyl(cage_center, root_position, cell, molname):
#     return Alkyl(cage_center, root_position, cell, molname, ["Ma", "Mb"])


# def _2_2_dimethylpropyl(cage_center, root_position, cell, molname):
#     """
#     2,2-dimethylpropyl group rooted at root_position toward cage_center.
#     """
#     return Alkyl(cage_center, root_position, cell, molname, ["Ma", ["Mb", "Mc", "Md", "Me"]])


# def _2_3_dimethylbutyl(cage_center, root_position, cell, molname):
#     """
#     put a butyl group rooted at root_position toward cage_center.
#     """
#     return Alkyl(
#         cage_center, root_position, cell, molname, [
#             "Ma", [
#                 "Mb", [
#                     "Mc", "Md", "Me"], "Mf"]])


# def _3_methylbutyl(cage_center, root_position, cell, molname):
#     """
#     put a butyl group rooted at root_position toward cage_center.
#     """
#     return Alkyl(cage_center, root_position, cell, molname, ["Ma", ["Mb", ["Mc", "Md", "Me"]]])


# def _3_3_dimethylbutyl(cage_center, root_position, cell, molname):
#     """
#     put a butyl group rooted at root_position toward cage_center.
#     """
#     return Alkyl(
#         cage_center, root_position, cell, molname, [
#             "Ma", [
#                 "Mb", [
#                     "Mc", "Md", "Me", "Mf"]]])


class GenIce():
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
    def __init__(self,
                 lat: Type[Lattice],
                 signature: str="",
                 density: float=0,
                 rep=(1, 1, 1),
                 cations:dict={},
                 anions:dict={},
                 spot_guests:dict={},
                 spot_groups:dict={},
                 asis:bool=False,
                 shift=(0., 0., 0.),
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
        # ================================================================
        # rotmatrices (analice)
        #
        try:
            self.rotmatrices = lat.rotmat
        except BaseException:
            logger.info("No rotmatrices in lattice")
            self.rotmatrices = None
        # ================================================================
        # cell: cell dimension
        #   see parse_cell for syntax.
        #
        self.cell = Cell(lat.cell)

        # ================================================================
        # waters: positions of water molecules
        #
        self.waters = put_in_array(lat.waters)
        logger.debug(f"Waters: {len(self.waters)}")
        self.waters = self.waters.reshape((self.waters.size // 3, 3))

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative", i.e. in fractional coordinates
        #
        if lat.coord == "absolute":
            self.waters = self.cell.abs2rel(self.waters)

        # shift of the origin
        self.waters = np.array(self.waters) + np.array(shift)
        # fractional coordinate between [0, 1)
        self.waters -= np.floor(self.waters)

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
            logger.info(f"Bond length (specified): {self.bondlen}")
        except AttributeError:
            logger.debug("  Estimating the bond threshold length...")
            # very rough estimate of the pair distances assuming that the particles distribute homogeneously.
            rc = (volume / nmol)**(1 / 3) * 1.5
            p = pl.pairs_iter(self.waters,
                              rc=rc,
                              cell=self.cell.mat,
                              distance=False)
            self.bondlen = 1.1 * \
                shortest_distance(self.waters, self.cell, pairs=p)
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
                    "Density is not specified. Assume the density from lattice.")
                dmin = shortest_distance(self.waters, self.cell)
                logger.info(
                    f"Closest pair distance: {dmin} (should be around 0.276 nm)")
                self.density = density0 / (0.276 / dmin)**3
                # self.density = density0
        else:
            self.density = density

        logger.info(f"Target Density: {self.density}")
        logger.info(f"Original Density: {density0}")

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density)**(1.0 / 3.0)
        self.cell.scale(ratio)

        if self.bondlen is not None:
            self.bondlen *= ratio
        logger.info(f"Bond length (scaled, nm): {self.bondlen}")

        # ================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        self.cagepos = None
        self.cagetype = None

        if "cages" in lat.__dict__:
            self.cagepos, self.cagetype = parse_cages(lat.cages)
            logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            # pre-parsed data
            self.cagepos, self.cagetype = np.array(lat.cagepos), lat.cagetype
        if self.cagepos is not None:
            self.cagepos = np.array(self.cagepos) + np.array(shift)
            self.cagepos -= np.floor(self.cagepos)

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

        if "dopeIonsToUnitCell" in lat.__dict__:
            self.dopeIonsToUnitCell = lat.dopeIonsToUnitCell
        else:
            self.dopeIonsToUnitCell = None
        self.dopants = set()

        # if asis, make pairs to be fixed.
        if self.asis and len(self.fixed) == 0:
            self.fixed = self.pairs

        # filled cages
        self.filled_cages = set()

        # groups info
        self.groups = defaultdict(dict)

        # groups for the semi-guest
        # experimental; there are many variation of semi-guest inclusion.
        # self.groups_placer = {"Bu-": butyl,
        #                       "Butyl-": butyl,
        #                       "Pentyl-": pentyl,
        #                       "Propyl-": propyl,
        #                       "2,2-dimethylpropyl-": _2_2_dimethylpropyl,
        #                       "2,3-dimethylbutyl-": _2_3_dimethylbutyl,
        #                       "3,3-dimethylbutyl-": _3_3_dimethylbutyl,
        #                       "3-methylbutyl-": _3_methylbutyl,
        #                       "Ethyl-": ethyl}

    def generate_ice(self,
                     formatter:Type[Format],
                     water:Union[Type[Molecule], None]=None,
                     guests={},
                     depol="strict",
                     noise=0.,
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
            logger.warn("generate_ice(water, formatter) is deprecated. "
                        "New syntax is: generate_ice(formatter, water=water).")

        def Stages():
            hooks = formatter.hooks()
            maxstage = max(0, *hooks.keys())

            if 0 in hooks:
                abort = hooks[0](self)
                if maxstage < 1 or abort:
                    return

            self.Stage1(noise,
                        assess_cages=assess_cages)

            if 1 in hooks:
                abort = hooks[1](self)
                if maxstage < 2 or abort:
                    return

            self.Stage2()

            # Count bonds
            num_hb_disorder = 0
            nfixed = 0
            for i, j, data in self.graph.edges(data=True):
                if self.graph[i][j]['fixed']:  # fixed pair
                    nfixed += 1
                else:
                    num_hb_disorder += 1
            logger.info(
                f"  Number of pre-oriented hydrogen bonds: {nfixed}")
            logger.info(
                f"  Number of unoriented hydrogen bonds: {num_hb_disorder}")
            logger.info("  Number of hydrogen bonds: {0} (regular num: {1})".format(
                nfixed + num_hb_disorder, len(self.reppositions) * 2))

            # test2==True means it is a z=4 graph.
            test2 = self.test_undirected_graph(self.graph)
            if not test2:
                logger.warn("Ice rule is not satisfied.")

            if 2 in hooks:
                abort = hooks[2](self)
                if maxstage < 3 or abort:
                    return

            # new_algorithm == fast algorithm
            new_algorithm = True
            # it makes the digraph obeying ice rule with zero net polarization
            # but it works only for a perfect 4-graph.
            if not test2 or self.asis or nfixed > 0 or depol != "strict":
                # The network is not 4-connected.
                new_algorithm = False

            if new_algorithm:
                # Fast track
                self.Stage3D()
            else:
                # Normal path; make it random, and then remove the defects.
                self.Stage3()

            if 3 in hooks:
                abort = hooks[3](self)
                if maxstage < 4 or abort:
                    return

            # spacegraph might be already set in Stage3D.
            if self.spacegraph is None:
                logger.debug(f"  graph? {self.spacegraph}")
                if num_hb_disorder == 0:
                    self.Stage4(depol="none")
                else:
                    self.Stage4(depol=depol)

            dipole = self.spacegraph.net_polarization()
            logger.info(f"Residual net polarization: {dipole[0]:.2f} {dipole[1]:.2f} {dipole[2]:.2f}")


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
    def Stage1(self,
               noise=0.,
               assess_cages=False):
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
        self.reppositions = replicate_positions(self.waters, self.rep)

        # This must be done before the replication of the cell.
        logger.info("  Number of water molecules: {0}".format(
            len(self.reppositions)))
        self.graph = self.prepare_random_graph(self.fixed)

        # scale the cell
        self.repcell = Cell(self.cell.mat)
        self.repcell.scale2(self.rep)

        if noise > 0.0:
            logger.info(f"  Add noise: {noise}.")
            perturb = np.random.normal(
                loc=0.0,
                scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                size=self.reppositions.shape)
            self.reppositions += self.repcell.abs2rel(perturb)

        if assess_cages:
            logger.info("  Assessing the cages...")
            self.cagepos, self.cagetype = cage.assess_cages( self.graph, self.waters )
            logger.info("  Done assessment.")

        if self.cagepos is not None and self.cagepos.shape[0] > 0:
            logger.info("  Hints:")
            self.repcagepos = replicate_positions(self.cagepos, self.rep)
            nrepcages = self.repcagepos.shape[0]
            self.repcagetype = [self.cagetype[i % len(self.cagetype)]
                                for i in range(nrepcages)]
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
        self.dopants = replicate_labeldict(
            self.dopants, len(self.waters), self.rep)
        self.groups = replicate_groups(
            self.groups, self.waters, self.cagepos, self.rep)

        # self.groups_info(self.groups)
        for root, cages in self.groups.items():
            self.filled_cages |= set(cages)

        # logger.info(("filled",self.filled_cages))
        # Replicate the graph
        self.graph = replicate_graph(self.graph, self.waters, self.rep)

        # Dope ions by options.
        if len(self.anions) > 0:
            logger.info(f"  Anionize: {self.anions}.")

            for site, name in self.anions.items():
                self.graph.anionize(site)
                self.dopants[site] = name

        if len(self.cations) > 0:
            logger.info(f"  Cationize: {self.cations}.")

            for site, name in self.cations.items():
                self.graph.cationize(site)
                self.dopants[site] = name

    @timeit
    @banner
    def Stage3(self):
        """
        Make a graph obeying the ice rule.

        Provided variables:
        graph: network obeying B-F rule.
        """

        logger = getLogger()

        if self.asis:
            logger.info("  Skip applying the ice rule by request.")
        else:
            self.graph.purge_ice_defects()

        self.spacegraph = None

    @timeit
    @banner
    def Stage4(self, depol="strict"):
        """
        Depolarize.

        Provided variables:
        spacegraph: depolarized network with node positions.
        """

        logger = getLogger()

        if self.asis:
            depol = "none"

        # self.spacegraph = dg.SpaceIceGraph(self.graph,
        #                                    coord=self.reppositions,
        #                                    immutables=self.graph.immutables)
        # dg.depolarize(self.spacegraph, self.repcell.mat, draw=None, depol=depol)
        digraph = dg.depolarize(self.graph,
                                coord=self.reppositions,
                                immutables=self.graph.immutables,
                                cell=self.repcell.mat,
                                depol=depol)
        # for debug
        # digraph = dg.depolarize(digraph,
        #                         coord=self.reppositions,
        #                         immutables=self.graph.immutables,
        #                         cell=self.repcell.mat,
        #                         depol=depol)
        self.spacegraph = dg.SpaceIceGraph(digraph,
                                           coord=self.reppositions,
                                           immutables=self.graph.immutables)

    @timeit
    @banner
    def Stage3D(self):
        """
        Tile the graph with directed cycles.
        """

        # Cに書きかえるなら、この下の3つをおきかえる。
        def cycle_edges(cycle):
            for i in range(len(cycle)):
                yield cycle[i - 1], cycle[i]

        @timeit
        @banner
        def spanningCycles(cycles):
            """
            Look up the traversal cycles.
            """
            dipoles = []
            spanning = []
            for j, cycle in enumerate(cycles):
                dipole = 0
                for a, b in cycle_edges(cycle):
                    # displacement vector
                    d = self.reppositions[b] - self.reppositions[a]
                    d -= np.floor(d + 0.5)
                    dipole += d
                if not np.allclose(dipole, 0):
                    # it is a cell-spanning cycle
                    dipoles.append(dipole)
                    spanning.append(j)
            dipoles = np.array(dipoles)
            return dipoles, spanning

        @timeit
        @banner
        def direct(dipoles, spanning):
            """
            Re-orient the cycles so as to minimize the net polarization.
            """
            bestm = 999999
            bestp = None
            dir = np.random.randint(2, size=len(dipoles)) * 2 - 1  # +1 or -1
            pol = dipoles.T @ dir
            pol2 = pol @ pol
            for i in range(len(dipoles) * 2):
                r = random.randint(0, len(dipoles) - 1)
                newpol = pol - 2 * dir[r] * dipoles[r]
                newpol2 = newpol @ newpol
                if newpol2 <= pol2:
                    dir[r] = -dir[r]
                    pol = newpol
                    pol2 = newpol2
                    if pol2 < 1e-6:
                        break
            logger.debug(f"  Depolarized to {pol} in {i} steps")
            return dir

        @timeit
        @banner
        def cycles2digraph(cycles):
            """
            Convert cycles to a digraph.
            """
            d = dg.IceGraph()
            for cycle in cycles:
                nx.add_cycle(d, cycle, fixed=False)
            return d

        logger = getLogger()

        pairs = np.array([(i, j)
                         for i, j in self.graph.edges()], dtype=np.int32)
        Nnode = len(self.reppositions)
        # cycles = [list(cycle) for cycle in tc.tile(pairs, Nnode, self.seed)]
        # Now uses the python version of tilecycles because it is fast enough.
        cycles = [list(cycle) for cycle in tc.tile(pairs, Nnode)]

        # evaluate the dipole of each cycle
        dipoles, spanning = spanningCycles(cycles)
        logger.debug(f"  {len(spanning)} spanning cycles.")

        # invert randomly to eliminate the net polarization.
        # Rarely, it cannot be depolarized.

        dir = direct(dipoles, spanning)
        pol = dipoles.T @ dir
        pol2 = pol @ pol

        # invert cycles
        for i, p in enumerate(dir):
            if p < 0:
                cycles[spanning[i]].reverse()

        import networkx as nx

        d = cycles2digraph(cycles)

        self.graph = d
        self.spacegraph = None
        if pol2 < 1e-6:
            # Skip Stage4
            logger.debug("  Depolarized in Stage3D.")
            self.spacegraph = dg.SpaceIceGraph(d,
                                               coord=self.reppositions,
                                               immutables=self.graph.immutables)

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
            self.reppositions, self.spacegraph, self.repcell)

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
        self.universe.append(arrange(self.reppositions,
                                     self.repcell,
                                     self.rotmatrices,
                                     water,
                                     immutables=set(self.dopants)))

    @timeit
    @banner
    def Stage7(self, guests):
        """
        Arrange guest atoms.

        Provided variables:
        atoms: atomic positions of all molecules.
        """

        logger = getLogger()

        if self.cagepos is not None:

            # the cages around the dopants.
            dopants_neighbors = self.dopants_info(
                self.dopants, self.reppositions, self.repcagepos, self.repcell)

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
                        r = random.randint(0, len(rooms) - 1)
                        room = rooms[r]

                        if resident[room] is None:
                            resident[room] = molec
                            molecules[molec].append(room)
                            movedin.append(room)
                            remain -= 1

            # Now ge got the address book of the molecules.
            if len(molecules):
                logger.info("  Summary of guest placements:")
                self.guests_info(self.cagetypes, molecules)

            if len(self.spot_groups) > 0:
                logger.info("  Summary of groups:")
                self.groups_info(self.groups)

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
                            cage_center, pos, self.repcell, molname, origin_atom=name))

            # molecular guests
            for molec, cages in molecules.items():
                guest_type, guest_options = plugin_option_parser(molec)
                logger.debug(f"Guest type: {guest_type}")
                gmol = safe_import("molecule", guest_type).Molecule(
                    **guest_options)

                try:
                    mdoc = gmol.__doc__.splitlines()
                except BaseException:
                    mdoc = []
                for line in mdoc:
                    logger.info("  " + line)
                cage_center = [self.repcagepos[i] for i in cages]
                cmat = [np.identity(3) for i in cages]
                self.universe.append(arrange(cage_center,
                                             self.repcell,
                                             cmat,
                                             gmol))

        # Assume the dopant is monatomic and replaces one water molecule
        atomset = defaultdict(set)
        for label, name in self.dopants.items():
            atomset[name].add(label)

        for name, labels in atomset.items():
            pos = [self.reppositions[i] for i in sorted(labels)]
            rot = [self.rotmatrices[i] for i in sorted(labels)]
            oneatom = one.Molecule(label=name)
            self.universe.append(arrange(pos,
                                         self.repcell,
                                         rot,
                                         oneatom))

    def prepare_random_graph(self, fixed):

        logger = getLogger()
        if self.pairs is None:
            logger.info("  Pairs are not given explicitly.")
            logger.info(
                "  Estimating the bonds according to the pair distances.")

            logger.debug(f"Bondlen: {self.bondlen}")
            # make bonded pairs according to the pair distance.
            # make before replicating them.
            grid = pl.determine_grid(self.cell.mat, self.bondlen)
            assert np.product(
                grid) > 0, "Too thin unit cell. Consider use of --rep option if the cell was made by cif2ice."
            self.pairs = pl.pairs_fine(
                self.waters, self.bondlen, self.cell.mat, grid, distance=False)

            # self.pairs = [v for v in zip(j0,j1)]
            # Check using a simpler algorithm.
            # Do not use it for normal debug because it is too slow
            if False:  # logger.level <= logging.DEBUG:
                pairs1 = self.pairs
                pairs2 = [v for v in pl.pairs_crude(
                    self.waters, self.bondlen, self.cell.mat, distance=False)]
                logger.debug(f"pairs1: {len(pairs1)}")
                logger.debug(f"pairs2: {len(pairs2)}")
                for pair in pairs1:
                    i, j = pair
                    assert (i, j) in pairs2 or (j, i) in pairs2
                for pair in pairs2:
                    i, j = pair
                    assert (i, j) in pairs1 or (j, i) in pairs1

        graph = dg.IceGraph()
        if fixed is not None:
            for i, j in fixed:
                graph.add_edge(i, j, fixed=True)

        # Fixed pairs are default.
        for pair in self.pairs:
            i, j = pair

            if graph.has_edge(i, j) or graph.has_edge(j, i):
                pass
            else:
                if random.randint(0, 1) == 0:
                    graph.add_edge(i, j, fixed=False)
                else:
                    graph.add_edge(j, i, fixed=False)

        logger.info(f"  Number of water nodes: {graph.number_of_nodes()}")

        return graph

    def test_undirected_graph(self, graph):
        # Test

        logger = getLogger()
        undir = graph.to_undirected()
        for node in range(undir.number_of_nodes()):
            if node not in undir:
                logger.debug(f"z=0 at {node}")
            else:
                z = len(list(undir.neighbors(node)))
                if z != 4:
                    logger.debug(f"z={z} at {node}")

        if graph.number_of_edges() != len(self.reppositions) * 2:
            logger.info(
                "Inconsistent number of HBs {0} for number of molecules {1}.".format(
                    graph.number_of_edges(), len(
                        self.reppositions)))
            return False

        return True

    def dopants_info(self, dopants=None, waters=None, cagepos=None, cell=None):
        logger = getLogger()
        if dopants is None:
            dopants = self.dopants

        if waters is None:
            waters = self.waters

        if cagepos is None:
            cagepos = self.cagepos

        if cell is None:
            cell = self.cell

        dopants_neighbors = neighbor_cages_of_dopants(
            dopants, waters, cagepos, cell)

        for dopant, cages in dopants_neighbors.items():
            logger.info(
                f"    Cages adjacent to dopant {dopant}: {cages}")

        return dopants_neighbors

    def groups_info(self, groups):
        logger = getLogger()
        for root, cages in groups.items():
            for cage, group in cages.items():
                logger.info(
                    f"    Group {group} of dopant {root} in cage {cage}")

    def guests_info(self, cagetypes, molecules):
        logger = getLogger()
        for cagetype, cageid in cagetypes.items():
            logger.info(f"    Guests in cage type {cagetype}:")

            for molec, cages in molecules.items():
                cages = set(cages)
                cages &= cageid

                if len(cages):
                    logger.info(
                        f"      {molec} * {len(cages)} @ {cages}")

    def add_group(self, cage, group, root):
        self.groups[root][cage] = group
        self.filled_cages.add(cage)

    def __del__(self):
        logger = getLogger()
        logger.info("Completed.")
