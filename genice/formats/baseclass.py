import logging
import random
import math
import itertools as it
import sys

import numpy as np
from collections import Iterable, defaultdict

from genice.importer import safe_import
from genice import digraph   as dg
from genice import pairlist as pl

def shortest_distance(coord, cell, pairs=None):
    dmin = 1e99
    if pairs is None:
        iter = it.combinations(coord,2)
    else:
        iter = [(coord[i],coord[j]) for i,j in pairs]
    for c1,c2 in iter:
        d = c1-c2
        d -= np.floor(d + 0.5)
        r = np.dot(d,cell)
        rr = np.dot(r,r)
        if rr < dmin:
            dmin = rr
    return dmin**0.5


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
        nei = graph.neighbors(node)
        oh1 = coord[nei[0]] - coord[node]
        oh1 -= np.floor( oh1 + 0.5 )
        oh1 = np.dot(oh1,cell)                #abs coord
        oh2 = coord[nei[1]] - coord[node]
        oh2 -= np.floor( oh2 + 0.5 )
        oh2 = np.dot(oh2,cell)                #abs coord
        #oh1 /= np.linalg.norm(oh1)
        #oh2 /= np.linalg.norm(oh2)
        #logger.debug("bond angle cos:{0}".format(np.dot(oh1,oh2)))
        y  = oh2 - oh1
        y /= np.linalg.norm(y)
        z = (oh1 + oh2)/2
        z /= np.linalg.norm(z)
        x = np.cross(y,z)
        rotmat = np.vstack([x,y,z])
        rotmatrices.append(rotmat)
    return rotmatrices



def arrange_atoms(coord, cell, rotmatrices, intra, labels, name):
    atoms = []
    if len(intra) == 0:
        return atoms
    for node in range(len(coord)):
        abscom = np.dot(coord[node],cell)  # relative to absolute
        rotated = np.dot(intra,rotmatrices[node])
        for i in range(len(labels)):
            atoms.append([i,name,labels[i],rotated[i,:]+abscom])
    return atoms



def replicate_graph(graph, positions, rep):
    repgraph = dg.IceGraph()
    nmol = positions.shape[0]
    for i,j in graph.edges_iter(data=False):
        vec = positions[j] - positions[i]    #positions in the unreplicated cell
        delta = np.floor( vec + 0.5 ).astype(int)
        for x in range(rep[0]):
            for y in range(rep[1]):
                for z in range(rep[2]):
                    xi = (x + delta[0] + rep[0]) % rep[0]
                    yi = (y + delta[1] + rep[1]) % rep[1]
                    zi = (z + delta[2] + rep[2]) % rep[2]
                    newi = i+nmol*(xi+rep[0]*(yi+rep[1]*zi))
                    newj = j+nmol*(x +rep[0]*(y +rep[1]*z ))
                    if graph[i][j]['fixed']: ##fixed pair
                        repgraph.add_edge(newi,newj,fixed=True)
                    else:
                        if 0 == random.randint(0,1):         #shuffle the bond directions
                            repgraph.add_edge(newi,newj,fixed=False)
                        else:
                            repgraph.add_edge(newj,newi,fixed=False)
    return repgraph



def flatten(item):
    """Yield items from any nested iterable; see REF."""
    if type(item) is str:
        yield item
    elif isinstance(item, Iterable):
        for x in item:
            yield from flatten(x)
    else:
        yield item
    

def replicate(positions, rep):
    repx = positions.copy()
    for m in range(1,rep[0]):
        v = np.array([m,0,0])
        repx = np.concatenate((repx,positions+v))
    repy = repx.copy()
    for m in range(1,rep[1]):
        v = np.array([0,m,0])
        repy = np.concatenate((repy,repx+v))
    repz = repy.copy()
    for m in range(1,rep[2]):
        v = np.array([0,0,m])
        repz = np.concatenate((repz,repy+v))
    return repz / rep
    


def parse_cages(cages):
    logger = logging.getLogger()
    cagetype = []
    cagepos = []
    if type(cages) is str:
        for line in cages.split("\n"):
            cols = line.split()
            if len(cols)>0:
                cagetype.append(cols[0])
                cagepos.append([float(x) for x in cols[1:4]])
        cagepos = np.array(cagepos)
    else:
        #Assume it is a list of list
        #flatten
        c = [x for x in flatten(cages)]
        while len(c):
            cagetype.append(c.pop(0))
            cagepos.append(c[:3])
            c = c[3:]
        cagepos = np.array(cagepos)
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
        beta = cell[3] * math.pi / 180.
        cell = np.array(((cell[0]*1.0, cell[1]*0.0, cell[2]*math.cos(beta)),
                             (cell[0]*0.0, cell[1]*1.0, cell[2]*0.0),
                             (cell[0]*0.0, cell[1]*0.0, cell[2]*math.sin(beta))))
        return cell.transpose()   #all the vector calculations are done in transposed manner.
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
        return np.reshape(cell, (3,3))
        #assert cell[0,1] == 0 and cell[0,2] == 0 and cell[1,2] == 0
    else:
        logger.error("unknown cell type: {0}".format(celltype))
        sys.exit(1)





class GenIce():
    def __init__(self, options):
        """
        setup with options
        load the lattice module
        """
        self.logger = logging.getLogger()
        #pick required data here
        lattice_type = options.Type[0]
        seed         = options.seed[0]
        self.rep     = options.rep
        self.density = options.dens[0]
        self.nodep   = options.nodep
        #Set random seeds
        random.seed(seed)
        np.random.seed(seed)
        #never access options from later on
        self.logger.info("Ice type: {0}".format(lattice_type))
        lat = safe_import("lattice", lattice_type)
        #Show the document of the module
        try:
            self.doc = lat.__doc__.splitlines()
        except:
            self.doc = []
        self.doc.append("")
        self.doc.append("Command line: {0}".format(" ".join(sys.argv)))
        for line in self.doc:
            self.logger.info("!!! {0}".format(line))
        #================================================================
        # waters: positions of water molecules
        #
        self.waters = load_numbers(lat.waters)
        self.logger.debug("Waters: {0}".format(len(self.waters)))
        self.waters = self.waters.reshape((self.waters.size//3,3))

        #================================================================
        # cell: cell dimension
        # celltype: symmetry of the cell
        #   see parse_cell for syntax.
        #
        self.cell = parse_cell(lat.cell, lat.celltype)
        
        #================================================================
        # coord: "relative" or "absolute"
        #   Inside genice, molecular positions are always treated as "relative"
        #
        if lat.coord == "absolute":
            self.waters = np.dot(self.waters,np.linalg.inv(self.cell),)
            self.waters = np.array(self.waters)

        #================================================================
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
                        i,j = [int(x) for x in columns]
                        self.pairs.append((i,j))
            elif type(lat.pairs) is list:
                self.pairs = lat.pairs
                #for pair in lat.pairs:
                #    self.pairs.append(pair)
        except AttributeError:
            self.logger.info("Graph is not defined.")

        #================================================================
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
        #Set density
        mass = 18  #water
        NB   = 6.022e23
        nmol = self.waters.shape[0]        #nmol in a unit cell
        volume = np.linalg.det(self.cell)  #volume of a unit cell in nm**3
        density0 = mass * nmol / (NB*volume*1e-21)
        if self.density <= 0:
            try:
                self.density = lat.density
            except AttributeError:
                self.logger.info("Density is not specified. Assume the density from lattice.")
                dmin = shortest_distance(self.waters, self.cell)
                self.logger.info("Closest pair distance: {0} (should be around 0.276 nm)".format(dmin))
                self.density = density0 / (0.276/dmin)**3
                #self.density = density0
        self.logger.info("Density: {0}".format(self.density))
        self.logger.info("Density0: {0}".format(density0))

        #scale the cell according to the (specified) density
        ratio = (density0 / self.density)**(1.0/3.0)
        self.cell *= ratio
        if self.bondlen is not None:
            self.bondlen   *= ratio
        self.logger.info("Bond length (scaled): {0}".format(self.bondlen))

        #================================================================
        # double_network: True or False
        #   This is a special option for ices VI and VII that have
        #   interpenetrating double network.
        #   GenIce's fast depolarization algorithm fails in some case.
        #
        try:
            self.double_network = lat.double_network
        except AttributeError:
            self.double_network = False

        #================================================================
        # cages: positions of the centers of cages
        #   In fractional coordinate.
        #
        try:
            cagepos, self.cagetype = parse_cages(lat.cages)
            self.cagepos = replicate(cagepos,self.rep)
        except AttributeError:
            self.cagepos = None
            self.cagetype = None
        self.nodep = options.nodep

        #================================================================
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
                        i,j = [int(x) for x in columns]
                        self.fixed.append((i,j))  #Is a tuple
            elif type(lat.fixed) is list:
                self.fixed = []
                for pair in lat.fixed:
                    self.fixed.append(tuple(pair[:2])) #Must be a tuple
        except AttributeError:
            self.fixed = []


    def stage1(self):
        """
        replicate water molecules to make a repeated cell
        """
        self.logger.info("Stage1: Replication.")
        self.reppositions = replicate(self.waters, self.rep)
        #This must be done before the replication of the cell.
        self.graph = self.prepare_random_graph(self.fixed)
        #scale the cell
        for d in range(3):
            self.cell[:,d] *= self.rep[d]
        self.logger.info("Stage1: end.")
    
    def stage2(self):
        """
        make a random graph and replicate.
        """
        self.logger.info("Stage2: Graph preparation.")
        self.graph = replicate_graph(self.graph, self.waters, self.rep)
        #test2==True means it is a z=4 graph.
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
        self.logger.info("Stage4: depolarization.")
        if self.nodep:
            self.logger.info("Skip depolarization by request.")
            self.yapresult = ""
            self.spacegraph = dg.SpaceIceGraph(self.graph,coord=self.reppositions)
        else:
            if self.double_network:
                if (self.rep[0] % 2 == 0) and (self.rep[1] % 2 == 0) and (self.rep[2] % 2 == 0):
                    pass
                else:
                    self.logger.error("In making the ice structure having the double network (e.g. ices 6 and 7), all the repetition numbers (--rep) must be even.")
                    sys.exit(1)
            self.spacegraph = dg.SpaceIceGraph(self.graph,coord=self.reppositions)
            draw = dg.YaplotDraw(self.reppositions, self.cell, data=self.spacegraph)
            self.yapresult  = dg.depolarize(self.spacegraph, self.cell, draw=draw)
        self.logger.info("Stage4: end.")

    def stage5(self):
        """
        orientations for rigid molecules.
        """
        self.logger.info("Stage5: Orientation.")
        #Add small noises to the molecular positions
        self.reppositions += np.dot(np.random.random(self.reppositions.shape)*0.01-0.005, np.linalg.inv(self.cell))
        #determine the orientations of the water molecules based on edge directions.
        self.rotmatrices = orientations(self.reppositions, self.spacegraph, self.cell)


        #Activate it.
        #logger.info("The network is not specified.  Water molecules will be orinented randomly.")
        #rotmatrices = [rigid.rand_rotation_matrix() for pos in positions]
        self.logger.info("Stage5: end.")
        

    def stage6(self, water_type):
        """
        arrange water atoms
        """
        self.logger.info("Stage6: Atomic positions of water.")
        #assert audit_name(water_type), "Dubious water name: {0}".format(water_type)
        #water = importlib.import_module("genice.molecules."+water_type)
        water = safe_import("molecule", water_type)
        self.atoms = arrange_atoms(self.reppositions, self.cell, self.rotmatrices, water.sites, water.labels, water.name)
        self.logger.info("Stage6: end.")
        
    def stage7(self, guests):
        """
        arrange guest atoms
        """
        self.logger.info("Stage7: Atomic positions of the guest.")
        if self.cagepos is not None:
            cagetypes = set(self.cagetype)
            self.logger.info("Cage types: {0}".format(cagetypes))
        if guests is not None and self.cagepos is not None:
            #Make the cage type to guest type correspondence
            guest_in_cagetype = dict()
            for arg in guests:
                key, value = arg[0].split("=")
                guest_in_cagetype[key] = value
            #replicate the cagetype array
            cagetype = np.array([self.cagetype[i%len(self.cagetype)] for i in range(self.cagepos.shape[0])])
            for ctype in cagetypes:
                #filter the cagepos
                cpos = self.cagepos[cagetype == ctype]
                #guest molecules are not rotated.
                cmat = np.array([np.identity(3) for i in range(cpos.shape[0])])
                #If the guest molecule type is given,
                if ctype in guest_in_cagetype:
                    gname = guest_in_cagetype[ctype]
                    #Always check before dynamic import
                    #assert audit_name(gname), "Dubious guest name: {0}".format(gname)
                    #gmol = importlib.import_module("genice.molecules."+gname)
                    gmol = safe_import("molecule", gname)
                    self.logger.info("{0} is in the cage type '{1}'".format(guest_in_cagetype[ctype], ctype))
                    self.atoms += arrange_atoms(cpos, self.cell, cmat, gmol.sites, gmol.labels, gmol.name)
        self.logger.info("Stage7: end.")


    def stage7B(self, guests):
        """
        arrange guest atoms
        put them in separate lists
        """
        self.logger.info("Stage7B: Atomic positions of the guest.")
        self.guestAtoms = defaultdict(list)
        self.nGuestAtoms = defaultdict(int)
        if self.cagepos is not None:
            cagetypes = set(self.cagetype)
            self.logger.info("Cage types: {0}".format(cagetypes))
        if guests is not None and self.cagepos is not None:
            #Make the cage type to guest type correspondence
            guest_in_cagetype = dict()
            for arg in guests:
                key, value = arg[0].split("=")
                guest_in_cagetype[key] = value
            #replicate the cagetype array
            cagetype = np.array([self.cagetype[i%len(self.cagetype)] for i in range(self.cagepos.shape[0])])
            for ctype in cagetypes:
                #filter the cagepos
                cpos = self.cagepos[cagetype == ctype]
                #guest molecules are not rotated.
                cmat = np.array([np.identity(3) for i in range(cpos.shape[0])])
                #If the guest molecule type is given,
                if ctype in guest_in_cagetype:
                    gname = guest_in_cagetype[ctype]
                    gmol = safe_import("molecule", gname)
                    self.logger.info("{0} is in the cage type '{1}'".format(guest_in_cagetype[ctype], ctype))
                    self.guestAtoms[gname] += arrange_atoms(cpos, self.cell, cmat, gmol.sites, gmol.labels, gmol.name)
                    self.nGuestAtoms[gname] += len(cpos)
        self.logger.info("Stage7B: end.")
                            
        
    def prepare_random_graph(self, fixed):
        if self.pairs is None:
            self.logger.info("  Pairs are not given explicitly.")
            self.logger.info("  Start estimating the bonds according to the pair distances.")
            #make bonded pairs according to the pair distance.
            #make before replicating them.
            grid = pl.determine_grid(self.cell, self.bondlen)
            self.pairs = [v for v in pl.pairlist_fine(self.waters, self.bondlen, self.cell, grid, distance=False)]
            #Check using a simpler algorithm.
            if self.logger.level <= logging.DEBUG:
                pairs2 = [v for v in pl.pairlist_crude(self.waters, self.bondlen, self.cell, distance=False)]
                self.logger.debug("pairs: {0}".format(len(self.pairs)))
                self.logger.debug("pairs2: {0}".format(len(pairs2)))
                for pair in self.pairs:
                    i,j = pair
                    assert (i,j) in pairs2 or (j,i) in pairs2
                for pair in pairs2:
                    i,j = pair
                    assert (i,j) in self.pairs or (j,i) in self.pairs

        graph = dg.IceGraph()
        for i,j in fixed:
            graph.add_edge(i,j,fixed=True)
        #Fixed pairs are default.
        for pair in self.pairs:
            i,j = pair
            if graph.has_edge(i,j) or graph.has_edge(j,i):
                pass
            else:
                if random.randint(0,1) == 0:
                    graph.add_edge(i,j,fixed=False)
                else:
                    graph.add_edge(j,i,fixed=False)
        return graph

    def test_undirected_graph(self, graph):
        #Test
        undir = graph.to_undirected()
        for node in range(undir.number_of_nodes()):
            if node not in undir:
                self.logger.debug("z=0 at {0}".format(node))
            else:
                z = len(undir.neighbors(node))
                if  z!= 4:
                    self.logger.debug("z={0} at {1}".format(z,node))
        if graph.number_of_edges() != len(self.reppositions)*2:
            self.logger.info("Inconsistent number of HBs {0} for number of molecules {1}.".format(graph.number_of_edges(),len(self.reppositions)))
            return False
        return True

    def __del__(self):
        self.logger.info("Completed.")
        
