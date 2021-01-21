from logging import getLogger
import argparse as ap
from collections import defaultdict
import random
import itertools as it
from textwrap import wrap, fill

import numpy as np
import pairlist as pl

from genice.importer import safe_import
from genice import digraph as dg
from genice import __version__
from genice.cell import rel_wrap, Cell
from genice.valueparsers import parse_cages, parse_pairs, put_in_array, flatten
import genice.plugins
from genice.decorators import timeit, banner

# for alkyl groups (Experimental)
from genice import alkyl


class SmartFormatter(ap.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            # return text[2:].splitlines()
            return [line for L in text[2:].splitlines() for line in wrap(L, width=55, drop_whitespace=False) ]
        # this is the RawTextHelpFormatter._split_lines
        return ap.HelpFormatter._split_lines(self, text, width)
    def _get_help_string(self, action):
        if callable(action.help):
            return action.help()
        return action.help


def descriptions(category):
    titles={ "lattice": {"system": "1. Lattice structures served with GenIce",
                         "extra":  "2. Lattice structures served by plugins",
                         "local":  "3. Lattice structures served locally",
                         "title":  "[Available lattice structures]"},
             "format": {"system": "1. Formatters served with GenIce",
                        "extra":  "2. Formatters served by plugins",
                        "local":  "3. Formatters served locally",
                        "title":  "[Available formatters]"},
             "loader": {"system": "1. File types served with GenIce",
                        "extra":  "2. File types served by plugins",
                        "local":  "3. File types served locally",
                        "title":  "[Available input file types]"},
             }
    mods = genice.plugins.scan(category)
    catalog = " \n \n{0}\n \n".format(titles[category]["title"])
    desc = mods["desc"]
    for group in ("system", "extra", "local"):
        desced = defaultdict(list)
        undesc = []
        for L in mods[group]:
            if L in desc:
                desced[desc[L]].append(L)
            else:
                undesc.append(L)
        for dd in desced:
            desced[dd] = ", ".join(desced[dd])
        catalog += "{0}\n \n".format(titles[category][group])
        table = ""
        for dd in sorted(desced, key=lambda x: desced[x]):
            table += "{0}\t{1}\n".format(desced[dd], dd)
        if table == "":
            table = "(None)\n"
        table = [fill(line, width=55, drop_whitespace=False, expand_tabs=True, tabsize=16, subsequent_indent=" "*16) for line in table.splitlines()]
        table = "\n".join(table)+"\n"
        undesc = " ".join(undesc)
        if undesc != "":
            undesc = "(Undocumented) " + undesc
        catalog += table + "----\n" + undesc + "\n \n \n"
    return catalog




# 遅延評価。descriptions()関数は重いので、必要なければ呼びたくない。
def help_type():
    return 'R|Crystal type (1c, 1h, etc. See https://github.com/vitroid/GenIce for available ice structures.)\n\nIf you want to analyze your own structures, please try analice tool.\n\n' + descriptions("lattice")

def help_format():
    return 'R|Specify the output file format. [gromacs]\n\n'+descriptions("format")

def getoptions():
    parser = ap.ArgumentParser(description='GenIce is a swiss army knife to generate hydrogen-disordered ice structures. (version {0})'.format(__version__), prog='genice', formatter_class=SmartFormatter)
    parser.add_argument('--version',
                        '-V',
                        action='version',
                        version='%(prog)s {0}'.format(__version__))
    parser.add_argument('--rep',
                        '-r',
                        nargs=3,
                        type=int,
                        dest='rep',
                        default=[1, 1, 1],
                        help='Repeat the unit cell along a, b, and c axes. [1,1,1]')
    parser.add_argument('--dens',
                        '-d',
                        type=float,
                        dest='dens',
                        default=-1,
                        help='Specify the ice density in g/cm3 (Guests are not included.)')
    parser.add_argument('--add_noise',
                        type=float,
                        dest='noise',
                        default=0.,
                        metavar='percent',
                        help='Add a Gauss noise with given width (SD) to the molecular positions of water. The value 1 corresponds to 1 percent of the molecular diameter of water.')
    parser.add_argument('--seed',
                        '-s',
                        type=int,
                        dest='seed',
                        default=1000,
                        help='Random seed [1000]')
    parser.add_argument('--format',
                        '-f',
                        dest='format',
                        default="gromacs",
                        metavar="name",
                        help=help_format)
    parser.add_argument('--water',
                        '-w',
                        dest='water',
                        default="tip3p",
                        metavar="model",
                        help='Specify water model. (tip3p, tip4p, etc.) [tip3p]')
    parser.add_argument('--guest',
                        '-g',
                        nargs=1,
                        dest='guests',
                        metavar="D=empty",
                        action="append",
                        help='Specify guest(s) in the cage type. (D=empty, T=co2*0.5+me*0.3, etc.)')
    parser.add_argument('--Guest',
                        '-G',
                        nargs=1,
                        dest='spot_guests',
                        metavar="13=me",
                        action="append",
                        help='Specify guest in the specific cage. (13=me, 32=co2, etc.)')
    parser.add_argument('--Group',
                        '-H',
                        nargs=1,
                        dest='groups',
                        metavar="13=bu-:0",
                        action="append",
                        help='Specify the group. (-H 13=bu-:0, etc.)')
    parser.add_argument('--anion',
                        '-a',
                        nargs=1,
                        dest='anions',
                        metavar="3=Cl",
                        action="append",
                        help='Specify a monatomic anion that replaces a water molecule. (3=Cl, 39=F, etc.)')
    parser.add_argument('--cation',
                        '-c',
                        nargs=1,
                        dest='cations',
                        metavar="3=Na",
                        action="append",
                        help='Specify a monatomic cation that replaces a water molecule. (3=Na, 39=NH4, etc.)')
    parser.add_argument('--visual',
                        dest='visual',
                        default="",
                        metavar="visual",
                        help='Specify the yaplot file to store the depolarization paths. [""]')
    parser.add_argument('--nodep',
                        action='store_true',
                        dest='nodep',
                        default=False,
                        help='No depolarization.')
    parser.add_argument('--asis',
                        action='store_true',
                        dest='asis',
                        default=False,
                        help='Assumes all given HB pairs to be fixed. No shuffle and no depolarization.')
    parser.add_argument('--debug',
                        '-D',
                        action='store_true',
                        dest='debug',
                        help='Output debugging info.')
    parser.add_argument('--quiet',
                        '-q',
                        action='store_true',
                        dest='quiet',
                        help='Do not output progress messages.')
    parser.add_argument('Type',
                        help=help_type)
    return parser.parse_args()


def assume_tetrahedral_vectors(v):
    """
    Assume missing vectors at a tetrahedral node.

    Given: known vectors.
    Returns: assumed vectors
    """

    assert len(v) > 0

    if len(v) == 3:

        return [-(v[0] + v[1] + v[2])]

    elif len(v) == 2:
        y = v[1] - v[0]
        y /= np.linalg.norm(y)
        z = v[1] + v[0]
        z /= np.linalg.norm(z)
        x = np.cross(y, z)
        v2 = (x * 8.0**0.5 - z) / 3.0
        v3 = (-x * 8.0**0.5 - z) / 3.0

        return [v2, v3]

    elif len(v) == 1:
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
    rotmatrices = []
    assert len(coord) == graph.number_of_nodes()  # just for a test of pure water

    for node in range(graph.number_of_nodes()):
        if node in graph.ignores:
            # for dopants; do not rotate
            rotmat = np.identity(3)
        else:
            vsucc = [cell.rel2abs(rel_wrap(coord[x] - coord[node])) for x in graph.successors(node)]

            if len(vsucc) < 2:  # TSL
                vpred = [cell.rel2abs(rel_wrap(coord[x] - coord[node])) for x in graph.predecessors(node)]
                vsucc = [x / np.linalg.norm(x) for x in vsucc]
                vpred = [x / np.linalg.norm(x) for x in vpred]

                if len(vpred) > 2:
                    vpred = vpred[:2]  # number of incoming bonds should be <= 2
                vcomp = assume_tetrahedral_vectors(vpred + vsucc)
                logger.debug("Node {0} vcomp {1} vsucc {2} vpred {3}".format(node, vcomp, vsucc, vpred))
                vsucc = (vsucc + vcomp)[:2]

            logger.debug("Node {0} vsucc {1}".format(node, vsucc))
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


def arrange_atoms(coord, cell, rotmatrices, intra, labels, name, ignores=set()):
    logger = getLogger()
    atoms = []

    if len(intra) == 0:
        return atoms

    for order, pos in enumerate(coord):
        if order in ignores:
            continue

        abscom = cell.rel2abs(pos)  # relative to absolute
        rotated = np.dot(intra, rotmatrices[order])

        for i in range(len(labels)):
            atoms.append([i, name, labels[i], rotated[i, :] + abscom, order])
        # logger.debug((np.linalg.norm(rotated[0] - rotated[1])))

    return atoms


def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99

    if pairs is None:
        iter = it.combinations(coord, 2)
    else:
        iter = [(coord[i], coord[j]) for i, j in pairs]

    for c1, c2 in iter:
        r = cell.rel2abs(rel_wrap(c1 - c2))
        rr = np.dot(r, r)

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
            sqdistance = np.dot(a, a)

            if sqdistance < 0.57**2:
                dnei[site].add(i)
                # logger.info((i,cagepos[i]))

    return dnei


# They should be separate plugins in the future.
def butyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", ["Mc", "Md"]]])


def pentyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", ["Mc", ["Md", "Me"]]]])


def propyl(cpos, root, cell, molname):
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", "Mc"]])


def ethyl(cpos, root, cell, molname):
    return Alkyl(cpos, root, cell, molname, ["Ma", "Mb"])


def _2_2_dimethylpropyl(cpos, root, cell, molname):
    """
    2,2-dimethylpropyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", "Mc", "Md", "Me"]])


def _2_3_dimethylbutyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", ["Mc", "Md", "Me"], "Mf"]])


def _3_methylbutyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", ["Mc", "Md", "Me"]]])


def _3_3_dimethylbutyl(cpos, root, cell, molname):
    """
    put a butyl group rooted at root toward cpos.
    """
    return Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", ["Mc", "Md", "Me", "Mf"]]])


def Alkyl(cpos, root, cell, molname, backbone):
    """
    put a normal-alkyl group rooted at root toward cpos.
    """
    logger = getLogger()
    # logger.info("  Put butyl at {0}".format(molname))
    v1abs = cell.rel2abs(rel_wrap(cpos - root))
    v1 = v1abs / np.linalg.norm(v1abs)

    origin = cell.rel2abs(root)
    CC = 0.154
    rawatoms = alkyl.alkyl(v1, v1abs * 1.5 / CC, backbone)

    atoms = []
    for i, atom in enumerate(rawatoms):
        atomname, pos = atom
        atompos = cell.abs_wrapf(pos * CC + origin)
        atoms.append([i, molname, atomname, atompos, 0])

    return atoms


class GenIce():
    def __init__(self,
                 lat,
                 argv,
                 density=0,
                 rep=(1, 1, 1),
                 cations=dict(),
                 anions=dict(),
                 spot_guests=dict(),
                 spot_groups=dict(),
                 asis=False,
                 ):

        self.logger = getLogger()
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

        self.doc.append("")
        self.doc.append("Command line: {0}".format(" ".join(argv)))

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

        if "cages" in lat.__dict__:
            self.cagepos, self.cagetype = parse_cages(lat.cages)
            self.logger.warn("Use of 'cages' in a lattice-plugin is deprecated.")
        elif "cagepos" in lat.__dict__:
            # pre-parsed data
            self.cagepos, self.cagetype = np.array(lat.cagepos), lat.cagetype

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
        self.groups_placer = {"Bu-": butyl,
                              "Butyl-": butyl,
                              "Pentyl-": pentyl,
                              "Propyl-": propyl,
                              "2,2-dimethylpropyl-": _2_2_dimethylpropyl,
                              "2,3-dimethylbutyl-": _2_3_dimethylbutyl,
                              "3,3-dimethylbutyl-": _3_3_dimethylbutyl,
                              "3-methylbutyl-": _3_methylbutyl,
                              "Ethyl-": ethyl}

    def generate_ice(self,
                     water_type,
                     guests,
                     formatter,
                     record_depolarization_path=None,
                     depolarize=True,
                     noise=0.):

        hooks = formatter.hooks
        arg   = formatter.arg
        maxstage = max(0, *hooks.keys())
        logger = getLogger()

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

        res = self.stage2()

        if 2 in hooks:
            abort = hooks[2](self)
            if maxstage < 3 or abort:
                return

        self.stage3()

        if 3 in hooks:
            abort = hooks[3](self)
            if maxstage < 4 or abort:
                return

        self.stage4(depolarize=depolarize,
                    record_depolarization_path=record_depolarization_path)

        if 4 in hooks:
            abort = hooks[4](self)
            if maxstage < 5 or abort:
                return

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

        self.stage7(guests)

        if 7 in hooks:
            hooks[7](self)


    @timeit
    @banner
    def stage1(self,
               noise=0.):
        """
        Replicate water molecules to make a repeated cell

        Provided variables:
        repposition: replicated molecular positions (CoM, relative)
        repcell:     replicated simulation cell shape matrix
        repcagetype: replicated cage types array
        repcagepos:  replicated cage positions (CoM, relative)
        cagetypes:   set of cage types
        """

        self.reppositions = replicate_positions(self.waters, self.rep)

        # This must be done before the replication of the cell.
        self.logger.info("  Number of water molecules: {0}".format(
            len(self.reppositions)))
        self.graph = self.prepare_random_graph(self.fixed)

        # scale the cell
        self.repcell = Cell(self.cell.mat)
        self.repcell.scale2(self.rep)

        if noise > 0.0:
            self.logger.info("  Add noise: {0}.".format(noise))
            perturb = np.random.normal(loc=0.0,
                                       scale=noise * 0.01 * 3.0 * 0.5,  # in percent, radius of water
                                       size=self.reppositions.shape)
            self.reppositions += self.repcell.abs2rel(perturb)

        if self.cagepos is not None:
            self.logger.info("  Hints:")
            self.repcagepos = replicate_positions(self.cagepos, self.rep)
            nrepcages = self.repcagepos.shape[0]
            self.repcagetype = [self.cagetype[i % len(self.cagetype)]
                                for i in range(nrepcages)]
            self.cagetypes = defaultdict(set)

            for i, typ in enumerate(self.repcagetype):
                self.cagetypes[typ].add(i)

            # INFO for cage types
            self.logger.info("    Cage types: {0}".format(list(self.cagetypes)))

            for typ, cages in self.cagetypes.items():
                self.logger.info("    Cage type {0}: {1}".format(typ, cages))
            # Up here move to stage 1.


    @timeit
    @banner
    def stage2(self):
        """
        Make a random graph and replicate.

        Provided variables:
        dopants:
        groups:  replicated positions of the chemical groups (CoM, relative)
        filled_cages:
        graph:   replicated network topology (bond orientation may be random)
        """


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
        self.logger.info("  Number of pre-oriented hydrogen bonds: {0}".format(nfixed))
        self.logger.info("  Number of unoriented hydrogen bonds: {0}".format(nrandom))
        self.logger.info("  Number of hydrogen bonds: {0} (regular num: {1})".format(nfixed + nrandom, len(self.reppositions) * 2))

        # test2==True means it is a z=4 graph.
        self.test2 = self.test_undirected_graph(self.graph)
        if not self.test2:
            self.logger.warn("Test2 failed.")

        return self.test2

    @timeit
    @banner
    def stage3(self):
        """
        Make a true ice graph.

        Provided variables:
        graph: network obeying B-F rule.
        """

        if self.asis:
            self.logger.info("  Skip applying the ice rule by request.")
        else:
            self.graph.purge_ice_defects()


    @timeit
    @banner
    def stage4(self, depolarize=True, record_depolarization_path=None):
        """
        Depolarize.

        Provided variables:
        spacegraph: depolarized network with node positions.
        yapresult:  Animation of the depolarization process in YaPlot format.
        """


        if not depolarize or self.asis:
            self.logger.info("  Skip depolarization by request. {0} {1}".format(depolarize, self.asis))
            self.yapresult = ""
            self.spacegraph = dg.SpaceIceGraph(self.graph,
                                               coord=self.reppositions,
                                               ignores=self.graph.ignores)
        else:
            if self.double_network:
                if (self.rep[0] % 2 == 0) and (self.rep[1] % 2 == 0) and (self.rep[2] % 2 == 0):
                    pass
                else:
                    self.logger.error("In making the ice structure having the double network (e.g. ices 6 and 7), all the repetition numbers (--rep) must be even.")
                    sys.exit(1)
            self.spacegraph = dg.SpaceIceGraph(self.graph,
                                               coord=self.reppositions,
                                               ignores=self.graph.ignores)
            if record_depolarization_path is not None:
                draw = dg.YaplotDraw(self.reppositions, self.repcell.mat, data=self.spacegraph)
                yapresult = dg.depolarize(self.spacegraph, self.repcell.mat, draw=draw)
                record_depolarization_path.write(yapresult)
            else:
                dg.depolarize(self.spacegraph, self.repcell.mat, draw=None)

    @timeit
    @banner
    def stage5(self):
        """
        Prepare orientations for rigid molecules.

        Provided variables:
        reppositions: molecular positions.
        rotmatrices:  rotation matrices for water molecules
        """

        # determine the orientations of the water molecules based on edge
        # directions.
        self.rotmatrices = orientations(
            self.reppositions, self.spacegraph, self.repcell)

        # Activate it.
        # logger.info("The network is not specified.  Water molecules will be orinented randomly.")
        # rotmatrices = [rigid.rand_rotation_matrix() for pos in positions]

    @timeit
    @banner
    def stage6(self, water_type):
        """
        Arrange water atoms and replacements

        Provided variables:
        atoms: atomic positions of water molecules. (absolute)
        """

        self.logger.info("Stage6: Atomic positions of water.")

        # assert audit_name(water_type), "Dubious water name: {0}".format(water_type)
        # water = importlib.import_module("genice.molecules."+water_type)
        water = safe_import("molecule", water_type)

        try:
            mdoc = water.__doc__.splitlines()
        except BaseException:
            mdoc = []

        for line in mdoc:
            self.logger.info("  "+line)

        self.atoms = arrange_atoms(self.reppositions,
                                   self.repcell,
                                   self.rotmatrices,
                                   water.sites,
                                   water.labels,
                                   water.name,
                                   ignores=set(self.dopants))


    @timeit
    @banner
    def stage7(self, guests):
        """
        Arrange guest atoms

        Provided variables:
        atoms: atomic positions of all molecules.
        """

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

            if guests is not None:

                # process the -g option
                for arg in guests:
                    self.logger.debug(arg[0])
                    cagetype, spec = arg[0].split("=")
                    assert cagetype in self.cagetypes, "Nonexistent cage type: {0}".format(cagetype)
                    resident = dict()
                    rooms = list(self.cagetypes[cagetype] - self.filled_cages)

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

            # Now ge got the address book of the molecules.
            if len(molecules):
                self.logger.info("  Summary of guest placements:")
                self.guests_info(self.cagetypes, molecules)

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
                self.atoms.append([0, molname, name, self.repcell.rel2abs(pos), 0])
                del self.dopants[root]  # processed.
                self.logger.debug((root, cages, name, molname, pos, rot))

                for cage, group in cages.items():
                    assert group in self.groups_placer
                    assert cage in dopants_neighbors[root]
                    cpos = self.repcagepos[cage]
                    self.atoms += self.groups_placer[group](cpos,
                                                            pos,
                                                            self.repcell,
                                                            molname)

            # molecular guests
            for molec, cages in molecules.items():
                gmol = safe_import("molecule", molec)

                try:
                    mdoc = gmol.__doc__.splitlines()
                except BaseException:
                    mdoc = []
                for line in mdoc:
                    logger.info("  "+line)
                cpos = [self.repcagepos[i] for i in cages]
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


    def prepare_random_graph(self, fixed):

        if self.pairs is None:
            self.logger.info("  Pairs are not given explicitly.")
            self.logger.info("  Estimating the bonds according to the pair distances.")

            self.logger.debug("Bondlen: {0}".format(self.bondlen))
            # make bonded pairs according to the pair distance.
            # make before replicating them.
            grid = pl.determine_grid(self.cell.mat, self.bondlen)
            assert np.product(grid) > 0, "Too thin unit cell. Consider use of --rep option if the cell was made by cif2ice."
            self.pairs = pl.pairs_fine(self.waters, self.bondlen, self.cell.mat, grid, distance=False)

            # self.pairs = [v for v in zip(j0,j1)]
            # Check using a simpler algorithm.
            # Do not use it for normal debug because it is too slow
            if False:  # self.logger.level <= logging.DEBUG:
                pairs1 = self.pairs
                pairs2 = [v for v in pl.pairs_crude(self.waters, self.bondlen, self.cell.mat, distance=False)]
                self.logger.debug("pairs1: {0}".format(len(pairs1)))
                self.logger.debug("pairs2: {0}".format(len(pairs2)))
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

        self.logger.info("  Number of water nodes: {0}".format(graph.number_of_nodes()))

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

        dopants_neighbors = neighbor_cages_of_dopants(dopants, waters, cagepos, cell)

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
