import sys
import logging
import random
import itertools as it
import logging
from math import sin, cos, pi
from collections import Iterable, defaultdict

import numpy as np

from genice.importer import safe_import
from genice import pairlist as pl
from genice import digraph as dg
from genice import rigid

def load_numbers(v):
    if type(v) is str:
        return np.fromstring(v, sep=" ")
    elif type(v) is list:
        return np.array(v)
    else:
        return v


def orientations(coord, graph, cell):
    """
    Does not work when two OHs are colinear
    """
    logger = logging.getLogger()
    rotmatrices = []
    for node in range(graph.number_of_nodes()):
        if node in graph.ignores:
            # for dopants; do not rotate
            rotmat = np.identity(3)
        else:
            nei = list(graph.neighbors(node))
            oh1 = coord[nei[0]] - coord[node]
            oh1 -= np.floor(oh1 + 0.5)
            oh1 = np.dot(oh1, cell)                # abs coord
            oh2 = coord[nei[1]] - coord[node]
            oh2 -= np.floor(oh2 + 0.5)
            oh2 = np.dot(oh2, cell)                # abs coord
            # oh1 /= np.linalg.norm(oh1)
            # oh2 /= np.linalg.norm(oh2)
            # logger.debug("bond angle cos:{0}".format(np.dot(oh1, oh2)))
            y = oh2 - oh1
            y /= np.linalg.norm(y)
            z = (oh1 + oh2) / 2
            z /= np.linalg.norm(z)
            x = np.cross(y, z)
            rotmat = np.vstack([x, y, z])
        rotmatrices.append(rotmat)
    return rotmatrices


def arrange_atoms(coord, cell, rotmatrices, intra, labels, name, ignores=set()):
    atoms = []
    if len(intra) == 0:
        return atoms
    for order, pos in enumerate(coord):
        if order in ignores:
            continue
        abscom = np.dot(pos, cell)  # relative to absolute
        rotated = np.dot(intra, rotmatrices[order])
        for i in range(len(labels)):
            atoms.append([i, name, labels[i], rotated[i, :] + abscom, order])
    return atoms


def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99
    if pairs is None:
        iter = it.combinations(coord, 2)
    else:
        iter = [(coord[i], coord[j]) for i, j in pairs]
    for c1, c2 in iter:
        d = c1 - c2
        d -= np.floor(d + 0.5)
        r = np.dot(d, cell)
        rr = np.dot(r, r)
        if rr < dmin:
            dmin = rr
    return dmin**0.5


def flatten(item):
    """Yield items from any nested iterable; see REF."""
    if type(item) is str:
        yield item
    elif isinstance(item, Iterable):
        for x in item:
            yield from flatten(x)
    else:
        yield item


def replicate_groups(groups, waters, cagepos, rep):
    """
    This is not that easy.
    """
    logger = logging.getLogger()
    # Storage for replicated groups
    newgroups = defaultdict(dict)
    for root, cages in groups.items():
        # Position of root (water) (fractional)
        root_pos = waters[root]
        for cage, group_name in cages.items():
            # Position of the cage (fractional)
            cage_pos = cagepos[cage]
            # Relative position of the cage
            delta = cage_pos - root_pos
            delta -= np.floor(delta + 0.5)
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


def replicate_graph(graph, positions, rep):
    repgraph = dg.IceGraph()
    nmol = positions.shape[0]
    for i, j in graph.edges(data=False):
        # positions in the unreplicated cell
        vec = positions[j] - positions[i]
        delta = np.floor(vec + 0.5).astype(int)
        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    xi = (x + delta[0] + rep[0]) % rep[0]
                    yi = (y + delta[1] + rep[1]) % rep[1]
                    zi = (z + delta[2] + rep[2]) % rep[2]
                    newi = i + nmol * (xi + rep[0] * (yi + rep[1] * zi))
                    newj = j + nmol * (x + rep[0] * (y + rep[1] * z))
                    if graph[i][j]['fixed']:  # fixed pair
                        repgraph.add_edge(newi, newj, fixed=True)
                    else:
                        # shuffle the bond directions
                        if 0 == random.randint(0, 1):
                            repgraph.add_edge(newi, newj, fixed=False)
                        else:
                            repgraph.add_edge(newj, newi, fixed=False)
    # replicate "ignores" == dopants list in the graph
    repgraph.ignores = replicate_labels(graph.ignores, nmol, rep)
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


def parse_cages(cages):
    logger = logging.getLogger()
    cagetype = []
    cagepos = []
    if type(cages) is str:
        for line in cages.split("\n"):
            cols = line.split()
            if len(cols) > 0:
                cagetype.append(cols[0])
                cagepos.append([float(x) for x in cols[1:4]])
        cagepos = np.array(cagepos)
    else:
        # Assume it is a list of list
        # flatten
        c = list(flatten(cages))
        while len(c):
            cagetype.append(c.pop(0))
            cagepos.append(c[:3])
            c = c[3:]
        cagepos = np.array(cagepos)
    cagepos[:] -= np.floor(cagepos[:])
    return cagepos, cagetype


def parse_cell(cell, celltype):
    logger = logging.getLogger()
    if celltype == "rect":
        if type(cell) is str:
            cell = np.fromstring(cell, sep=" ")
        elif type(cell) is list:
            cell = np.array(cell)
        logger.debug("parse_cell 1: {0}".format(cell))
        return np.diag(cell)
    elif celltype == "monoclinic":
        if type(cell) is str:
            cell = np.fromstring(cell, sep=" ")
        elif type(cell) is list:
            cell = np.array(cell)
        beta = cell[3] * pi / 180.
        cell = np.array(((cell[0] * 1.0, cell[1] * 0.0, cell[2] * cos(beta)),
                         (cell[0] * 0.0, cell[1] * 1.0, cell[2] * 0.0),
                         (cell[0] * 0.0, cell[1] * 0.0, cell[2] * sin(beta))))
        # all the vector calculations are done in transposed manner.
        return cell.transpose()
    elif celltype == "triclinic":
        """
        Put the vectors like following:
        cell = "ax 0 0 bx by 0 cx cy cz"
        when you define a unit cell in Lattice/

        """
        if type(cell) is str:
            cell = np.fromstring(cell, sep=" ")
            logger.debug(cell)
        elif type(cell) is list:
            cell = np.array(cell)
        return np.reshape(cell, (3, 3))
        # assert cell[0, 1] == 0 and cell[0, 2] == 0 and cell[1, 2] == 0
    else:
        logger.error("unknown cell type: {0}".format(celltype))
        sys.exit(1)


def neighbor_cages_of_dopants(dopants, waters, cagepos, cell):
    """
    Just shows the environments of the dopants
    """
    #logger = logging.getLogger()
    dnei = defaultdict(set)
    for site, name in dopants.items():
        org = waters[site]
        for i, pos in enumerate(cagepos):
            #Displacement (relative)
            d = pos - org
            d -= np.floor(d + 0.5)
            #Displacement (absolute)
            a = np.dot(d, cell)
            sqdistance = np.dot(a, a)
            if sqdistance < 0.57**2:
                dnei[site].add(i)
                # logger.info((i,cagepos[i]))
    return dnei


def butyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    logger = logging.getLogger()
    # logger.info("  Put butyl at {0}".format(molname))
    v1 = cpos - root
    v1 -= np.floor(v1 + 0.5)
    v1 = np.dot(v1, cell)
    v1 /= np.linalg.norm(v1)
    v2 = np.random.random(3)
    v12 = np.dot(v1, v2)
    v2 -= v1 * v12
    v2 /= np.linalg.norm(v2)  # a random unit vector perpendicular to v1
    logger.debug("  {0} {1} {2}".format(
        np.dot(v1, v1), np.dot(v2, v2), np.dot(v1, v2)))
    #v3 = np.cross(v1,v2)
    origin = np.dot(root, cell)
    CC = 0.154
    c = cos(109.5 / 2 * pi / 180)
    s = (1.0 - c**2)**0.5
    atoms = []
    for i, atom in enumerate(["Ma", "Mb", "Mc", "Md"]):
        x = CC * s * (i + 1)
        y = ((i + 1) % 2) * CC * c
        atompos = x * v1 + y * v2 + origin
        atoms.append([i, molname, atom, atompos, 0])
    return atoms


class Lattice():
    def __init__(self,
                 lattice_type=None,
                 density=0,
                 rep=(1, 1, 1),
                 depolarize=True,
                 cations=dict(),
                 anions=dict(),
                 spot_guests=dict(),
                 spot_groups=dict(),
                 ):
        self.logger      = logging.getLogger()
        self.lattice_type = lattice_type
        self.rep         = rep
        self.depolarize  = depolarize
        self.cations     = cations
        self.anions      = anions
        self.spot_guests = spot_guests
        self.spot_groups = spot_groups
        lat = safe_import("lattice", lattice_type)
        # Show the document of the module
        try:
            self.doc = lat.__doc__.splitlines()
        except:
            self.doc = []
        self.doc.append("")
        self.doc.append("Command line: {0}".format(" ".join(sys.argv)))
        for line in self.doc:
            self.logger.info("!!! {0}".format(line))
        # ================================================================
        # waters: positions of water molecules
        #
        self.waters = load_numbers(lat.waters)
        self.logger.debug("Waters: {0}".format(len(self.waters)))
        self.waters = self.waters.reshape((self.waters.size // 3, 3))

        # ================================================================
        # cell: cell dimension
        # celltype: symmetry of the cell
        #   see parse_cell for syntax.
        #
        self.cell = parse_cell(lat.cell, lat.celltype)

        # ================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative"
        #
        if lat.coord == "absolute":
            self.waters = np.dot(self.waters, np.linalg.inv(self.cell), )
        self.waters = np.array([w - np.floor(w) for w in self.waters])

        # ================================================================
        # pairs: specify the pairs of molecules that are connected.
        #   Bond orientation will be shuffled later
        #   unless it is "fixed".
        #
        self.pairs = None
        try:
            if type(lat.pairs) is str:
                lines = lat.pairs.split("\n")
                self.pairs = []
                for line in lines:
                    columns = line.split()
                    if len(columns) == 2:
                        i, j = [int(x) for x in columns]
                        self.pairs.append((i, j))
            elif type(lat.pairs) is list:
                self.pairs = lat.pairs
                # for pair in lat.pairs:
                #    self.pairs.append(pair)
        except AttributeError:
            self.logger.info("Graph is not defined.")

        # ================================================================
        # bondlen: specify the bond length threshold.
        #   This is used when "pairs" are not specified.
        #   It is applied to the original positions of molecules (before density setting).
        #
        self.bondlen = None
        try:
            self.bondlen = lat.bondlen
            self.logger.info("Bond length: {0}".format(self.bondlen))
        except AttributeError:
            self.bondlen = 1.1 * shortest_distance(self.waters, self.cell)
            self.logger.info("Bond length (assumed): {0}".format(self.bondlen))
        # Set density
        mass = 18  # water
        NB = 6.022e23
        nmol = self.waters.shape[0]        # nmol in a unit cell
        volume = np.linalg.det(self.cell)  # volume of a unit cell in nm**3
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
        self.logger.info("Density: {0}".format(self.density))
        self.logger.info("Density0: {0}".format(density0))

        # scale the cell according to the (specified) density
        ratio = (density0 / self.density)**(1.0 / 3.0)
        self.cell *= ratio
        if self.bondlen is not None:
            self.bondlen *= ratio
        self.logger.info("Bond length (scaled): {0}".format(self.bondlen))

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
        if "cages" in lat.__dict__:
            self.cagepos, self.cagetype = parse_cages(lat.cages)

        # ================================================================
        # fixed: specify the bonds whose directions are fixed.
        #   you can specify them in pairs at a time.
        #   You can also leave it undefined.
        #
        try:
            if type(lat.fixed) is str:
                lines = lat.fixed.split("\n")
                self.fixed = []
                for line in lines:
                    columns = line.split()
                    if len(columns) == 2:
                        i, j = [int(x) for x in columns]
                        self.fixed.append((i, j))  # Is a tuple
            elif type(lat.fixed) is list:
                self.fixed = []
                for pair in lat.fixed:
                    self.fixed.append(tuple(pair[:2]))  # Must be a tuple
        except AttributeError:
            self.fixed = []
        if "dopeIonsToUnitCell" in lat.__dict__:
            self.dopeIonsToUnitCell = lat.dopeIonsToUnitCell
        else:
            self.dopeIonsToUnitCell = None
        self.dopants = set()

        # filled cages
        self.filled_cages = set()
        # groups info
        self.groups = defaultdict(dict)

        # groups for the semi-guest
        self.groups_placer = dict()
        self.groups_placer["Bu-"] = butyl  # is a function


    def format(self, water_type, guests, formatter):
        if 0 in formatter.hooks:
            formatter.hooks[0](self)
        if max(0,*formatter.hooks.keys()) < 1:
            return
        self.stage1()
        if 1 in formatter.hooks:
            formatter.hooks[1](self)
        if max(0,*formatter.hooks.keys()) < 2:
            return
        res = self.stage2()
        if 2 in formatter.hooks:
            formatter.hooks[2](self)
        if max(0,*formatter.hooks.keys()) < 3:
            return
        if not res:
            self.rotmatrices = [rigid.rand_rotation_matrix() for pos in self.reppositions]
        else:
            self.stage3()
            if 3 in formatter.hooks:
                formatter.hooks[3](self)
            if max(0,*formatter.hooks.keys()) < 4:
                return
            self.stage4()
            if 4 in formatter.hooks:
                formatter.hooks[4](self)
            if max(0,*formatter.hooks.keys()) < 5:
                return
            self.stage5()
            if 5 in formatter.hooks:
                formatter.hooks[5](self)
        if max(0,*formatter.hooks.keys()) < 6:
            return
        self.stage6(water_type)
        if 6 in formatter.hooks:
            formatter.hooks[6](self)
        if max(0,*formatter.hooks.keys()) < 7:
            return
        self.stage7(guests)
        if 7 in formatter.hooks:
            formatter.hooks[7](self)

        
    def stage1(self):
        """
        replicate water molecules to make a repeated cell
        """
        self.logger.info("Stage1: Replication.")
        self.reppositions = replicate_positions(self.waters, self.rep)
        # This must be done before the replication of the cell.
        self.logger.info("  Number of water molecules: {0}".format(
            len(self.reppositions)))
        self.graph = self.prepare_random_graph(self.fixed)
        # scale the cell
        self.repcell = np.zeros_like(self.cell)
        for d in range(3):
            self.repcell[d, :] = self.cell[d, :] * self.rep[d]
        self.logger.info("Stage1: end.")

    def stage2(self):
        """
        make a random graph and replicate.
        """
        self.logger.info("Stage2: Graph preparation.")
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
        # self.logger.info(("filled",self.filled_cages))
        # Replicate the graph
        self.graph = replicate_graph(self.graph, self.waters, self.rep)
        # Dope ions by options.
        if len(self.anions) > 0:
            self.logger.info("  Anionize: {0}.".format(self.anions))
            for site, name in self.anions.items():
                self.graph.anionize(site)
                self.dopants[site] = name
        if len(self.cations) > 0:
            self.logger.info("  Cationize: {0}.".format(self.cations))
            for site, name in self.cations.items():
                self.graph.cationize(site)
                self.dopants[site] = name
        # Count bonds
        nrandom = 0
        nfixed = 0
        for i, j, data in self.graph.edges(data=True):
            if self.graph[i][j]['fixed']:  # fixed pair
                nfixed += 1
            else:
                nrandom += 1
        self.logger.info(
            "  Number of pre-oriented hydrogen bonds: {0}".format(nfixed))
        self.logger.info(
            "  Number of unoriented hydrogen bonds: {0}".format(nrandom))
        self.logger.info("  Number of hydrogen bonds: {0} (regular num: {1})".format(
            nfixed + nrandom, len(self.reppositions) * 2))
        # test2==True means it is a z=4 graph.
        self.test2 = self.test_undirected_graph(self.graph)
        if not self.test2:
            self.logger.info("Test2 failed.")

        self.logger.info("Stage2: end.")
        return self.test2

    def stage3(self):
        """
        make a true ice graph.
        """
        self.logger.info("Stage3: Bernal-Fowler rule.")
        self.graph.purge_ice_defects()
        self.logger.info("Stage3: end.")

    def stage4(self):
        """
        depolarize.
        """
        self.logger.info("Stage4: Depolarization.")
        if not self.depolarize:
            self.logger.info("  Skip depolarization by request.")
            self.yapresult = ""
            self.spacegraph = dg.SpaceIceGraph(self.graph,
                                               coord=self.reppositions,
                                               ignores=self.graph.ignores)
        else:
            if self.double_network:
                if (self.rep[0] % 2 == 0) and (self.rep[1] % 2 == 0) and (self.rep[2] % 2 == 0):
                    pass
                else:
                    self.logger.error(
                        "In making the ice structure having the double network (e.g. ices 6 and 7), all the repetition numbers (--rep) must be even.")
                    sys.exit(1)
            self.spacegraph = dg.SpaceIceGraph(self.graph,
                                               coord=self.reppositions,
                                               ignores=self.graph.ignores)
            draw = dg.YaplotDraw(
                self.reppositions, self.repcell, data=self.spacegraph)
            self.yapresult = dg.depolarize(
                self.spacegraph, self.repcell, draw=draw)
        self.logger.info("Stage4: end.")

    def stage5(self):
        """
        orientations for rigid molecules.
        """
        self.logger.info("Stage5: Orientation.")
        # Add small noises to the molecular positions
        self.reppositions += np.dot(np.random.random(self.reppositions.shape)
                                    * 0.01 - 0.005, np.linalg.inv(self.repcell))
        # determine the orientations of the water molecules based on edge
        # directions.
        self.rotmatrices = orientations(
            self.reppositions, self.spacegraph, self.repcell)

        # Activate it.
        # logger.info("The network is not specified.  Water molecules will be orinented randomly.")
        # rotmatrices = [rigid.rand_rotation_matrix() for pos in positions]
        self.logger.info("Stage5: end.")

    def stage6(self, water_type):
        """
        arrange water atoms and replacements
        """
        self.logger.info("Stage6: Atomic positions of water.")
        # assert audit_name(water_type), "Dubious water name: {0}".format(water_type)
        # water = importlib.import_module("genice.molecules."+water_type)
        water = safe_import("molecule", water_type)
        self.atoms = arrange_atoms(self.reppositions,
                                   self.repcell,
                                   self.rotmatrices,
                                   water.sites,
                                   water.labels,
                                   water.name,
                                   ignores=set(self.dopants))
        self.logger.info("Stage6: end.")

    def stage7(self, guests):
        """
        arrange guest atoms
        """
        self.logger.info("Stage7: Atomic positions of the guest.")
        if self.cagepos is not None:
            self.logger.info("  Hints:")
            repcagepos = replicate_positions(self.cagepos, self.rep)
            repcagetype = [self.cagetype[i % len(self.cagetype)]
                           for i in range(repcagepos.shape[0])]
            cagetypes = defaultdict(set)
            for i, typ in enumerate(repcagetype):
                cagetypes[typ].add(i)
            # INFO for cage types
            self.logger.info("    Cage types: {0}".format(list(cagetypes)))
            for typ, cages in cagetypes.items():
                self.logger.info("    Cage type {0}: {1}".format(typ, cages))
            # the cages around the dopants.
            dopants_neighbors = self.dopants_info(
                self.dopants, self.reppositions, repcagepos, self.repcell)
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
            if guests is not None:
                # process the -g option
                for arg in guests:
                    cagetype, spec = arg[0].split("=")
                    assert cagetype in cagetypes, "Nonexistent cage type: {0}".format(
                        cagetype)
                    resident = dict()
                    rooms = list(cagetypes[cagetype] - self.filled_cages)
                    for room in rooms:
                        resident[room] = None
                    # spec contains a formula consisting of "+" and "*"
                    contents = spec.split("+")
                    vacant = len(rooms)
                    for content in contents:
                        if "*" in content:
                            molec, frac = content.split("*")
                            frac = float(frac)
                        else:
                            molec = content
                            frac = 1.0
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
                        #self.logger.info(
                        #    "    {0} * {1} @ {2}".format(molec, nmolec, movedin))
            # Now ge got the address book of the molecules.
            if len(molecules):
                self.logger.info("  Summary of guest placements:")
                self.guests_info(cagetypes, molecules)
            if len(self.spot_groups) > 0:
                self.logger.info("  Summary of groups:")
                self.groups_info(self.groups)
            # semi-guests
            for root, cages in self.groups.items():
                assert root in self.dopants
                name = self.dopants[root]
                molname = "G{0}".format(root)
                pos = self.reppositions[root]
                rot = self.rotmatrices[root]
                self.atoms.append([0, molname, name, np.dot(pos, self.repcell), 0])
                del self.dopants[root]  # processed.
                self.logger.debug((root,cages,name,molname,pos,rot))
                for cage, group in cages.items():
                    assert group in self.groups_placer
                    assert cage in dopants_neighbors[root]
                    cpos = repcagepos[cage]
                    self.atoms += self.groups_placer[group](cpos,
                                                            self.reppositions[root],
                                                            self.repcell,
                                                            molname)
            # molecular guests
            for molec, cages in molecules.items():
                gmol = safe_import("molecule", molec)
                cpos = [repcagepos[i] for i in cages]
                cmat = [np.identity(3) for i in cages]
                self.atoms += arrange_atoms(cpos, self.repcell,
                                            cmat, gmol.sites, gmol.labels, gmol.name)
        # Assume the dopant is monatomic and replaces one water molecule
        atomset = defaultdict(set)
        for label, name in self.dopants.items():
            atomset[name].add(label)
        for name, labels in atomset.items():
            pos = [self.reppositions[i] for i in sorted(labels)]
            rot = [self.rotmatrices[i] for i in sorted(labels)]
            self.atoms += arrange_atoms(pos,
                                        self.repcell,
                                        rot,
                                        [[0., 0., 0.], ],
                                        [name],
                                        name)
        self.logger.info("Stage7: end.")


    def prepare_random_graph(self, fixed):
        if self.pairs is None:
            self.logger.info("  Pairs are not given explicitly.")
            self.logger.info(
                "  Start estimating the bonds according to the pair distances.")
            # make bonded pairs according to the pair distance.
            # make before replicating them.
            grid = pl.determine_grid(self.cell, self.bondlen)
            assert np.product(grid) > 0, "Too thin unit cell. Consider use of --rep option if the cell was made by cif2ice."
            self.pairs = [v for v in pl.pairlist_fine(
                self.waters, self.bondlen, self.cell, grid, distance=False)]
            # Check using a simpler algorithm.
            if self.logger.level <= logging.DEBUG:
                pairs2 = [v for v in pl.pairlist_crude(
                    self.waters, self.bondlen, self.cell, distance=False)]
                self.logger.debug("pairs: {0}".format(len(self.pairs)))
                self.logger.debug("pairs2: {0}".format(len(pairs2)))
                for pair in self.pairs:
                    i, j = pair
                    assert (i, j) in pairs2 or (j, i) in pairs2
                for pair in pairs2:
                    i, j = pair
                    assert (i, j) in self.pairs or (j, i) in self.pairs

        graph = dg.IceGraph()
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
        return graph

    def test_undirected_graph(self, graph):
        # Test
        undir = graph.to_undirected()
        for node in range(undir.number_of_nodes()):
            if node not in undir:
                self.logger.debug("z=0 at {0}".format(node))
            else:
                z = len(list(undir.neighbors(node)))
                if z != 4:
                    self.logger.debug("z={0} at {1}".format(z, node))
        if graph.number_of_edges() != len(self.reppositions) * 2:
            self.logger.info("Inconsistent number of HBs {0} for number of molecules {1}.".format(
                graph.number_of_edges(), len(self.reppositions)))
            return False
        return True

    def dopants_info(self, dopants=None, waters=None, cagepos=None, cell=None):
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
            self.logger.info(
                "    Cages adjacent to dopant {0}: {1}".format(dopant, cages))
        return dopants_neighbors

    def groups_info(self, groups):
        for root, cages in groups.items():
            for cage, group in cages.items():
                self.logger.info(
                    "    Group {0} of dopant {1} in cage {2}".format(group, root, cage))

    def guests_info(self, cagetypes, molecules):
        for cagetype, cageid in cagetypes.items():
            self.logger.info("    Guests in cage type {0}:".format(cagetype))
            for molec, cages in molecules.items():
                cages = set(cages)
                cages &= cageid
                if len(cages):
                    self.logger.info("      {0} * {1} @ {2}".format(molec, len(cages), cages))
        
    def add_group(self, cage, group, root):
        self.groups[root][cage] = group
        self.filled_cages.add(cage)

    def __del__(self):
        self.logger.info("Completed.")
