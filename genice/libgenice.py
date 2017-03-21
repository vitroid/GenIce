#!/usr/bin/env python
# -*- python -*-

from __future__ import print_function
import sys
import numpy     as np
import random
import math
import logging
import re
import importlib
import os
from genice import pairlist  as pl
from genice import digraph   as dg
import networkx as nx
from collections import defaultdict
from collections import Iterable

def audit_name(name):
    """
    Audit the mol name to avoid the access to external files
    """
    return re.match('^[A-Za-z0-9-_]+$', name) is not None



def safe_import(category, name):
    assert category in ("lattice", "molecule")
    assert audit_name(name), "Dubious {0} name: {1}".format(category, name)
    module = None
    if category == "lattice":
        try:
            module     = importlib.import_module("lattices."+name) #at ~/.genice
        except ImportError as e:
            pass
        if module is None:
            module     = importlib.import_module("genice.lattices."+name)
    else:
        try:
            module     = importlib.import_module("molecules."+name) #at ~/.genice
        except ImportError:
            pass
        if module is None:
            module     = importlib.import_module("genice.molecules."+name)
    return module


def usage(parser):
    parser.print_help()
    sys.exit(1)



def orientations(coord, graph, cell):
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



def replicate_graph(graph, positions, rep):
    repgraph = dg.IceGraph()
    nmol = positions.shape[0]
    for i,j,k in graph.edges_iter(data=True):
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
                    if 0 == random.randint(0,1):         #shuffle the bond directions
                        repgraph.add_edge(newi,newj)
                    else:
                        repgraph.add_edge(newj,newi)
    return repgraph



def load_numbers(v):
    if type(v) is str:
        return np.fromstring(v, sep=" ")
    elif type(v) is list:
        return np.array(v)
    

def flatten(item):
    """Yield items from any nested iterable; see REF."""
    if type(item) is str:
        yield item
    elif isinstance(item, Iterable):
        for x in item:
            yield from flatten(x)
    else:
        yield item
    
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


#0th stage: determine the molecular positions
#1st stage: undirected graph (connectivity)
#2nd stage: directed graph (obeying ice rule)
#3rd stage: depolarized graph


def generate_ice(lattice_type, density=-1, seed=1000, rep=(1,1,1), stage=3):
    logger = logging.getLogger()
    logger.info("Ice type: {0}".format(lattice_type))
    lat = safe_import("lattice", lattice_type)
    #Show the document of the module
    try:
        for line in lat.__doc__.splitlines():
            logger.info("!!! {0}".format(line))
    except:
        pass
    lat.waters = load_numbers(lat.waters)
    logger.debug("Waters: {0}".format(len(lat.waters)))
    lat.waters = lat.waters.reshape((lat.waters.size//3,3))
    #prepare cell transformation matrix
    lat.cell = parse_cell(lat.cell, lat.celltype)
    #express molecular positions in the coordinate relative to the cell
    
    if lat.coord == "absolute":
        lat.waters = np.dot(lat.waters,np.linalg.inv(lat.cell),)
        lat.waters = np.array(lat.waters)
    random.seed(seed)
    np.random.seed(seed)

    #Prearranged network topology information (if available)
    pairs = None
    bondlen = None
    try:
        if type(lat.pairs) is str:
            lines = lat.pairs.split("\n")
            pairs = set()
            for line in lines:
                columns = line.split()
                if len(columns) == 2:
                    i,j = [int(x) for x in columns]
                    pairs.add(frozenset([i,j]))
        elif type(lat.pairs) is list:
            pairs = set()
            for pair in pairs:
                pairs.add(frozenset(pair))
    except AttributeError:
        logger.info("Graph is not defined.")

    #Bond length threshold
    try:
        bondlen = lat.bondlen
    except AttributeError:
        logger.info("Bond length is not defined.")


    #set density
    if density < 0:
        density = lat.density

    #scale the cell according to the (specified) density
    mass = 18  #water
    NB   = 6.022e23
    nmol = lat.waters.shape[0]        #nmol in a unit cell
    volume = np.linalg.det(lat.cell)  #volume of a unit cell in nm**3
    d    = mass * nmol / (NB*volume*1e-21)
    ratio = (d / density)**(1.0/3.0)
    lat.cell *= ratio
    if bondlen is not None:
        bondlen   *= ratio

    if True:
        logger.info("Start placing the bonds.")
        if pairs is None:
            logger.info("  Pairs are not given explicitly.")
            logger.info("  Start estimating the bonds according to the pair distances.")
            #make bonded pairs according to the pair distance.
            #make before replicating them.
            grid = pl.determine_grid(lat.cell, bondlen)
            pairs = [v for v in pl.pairlist_fine(lat.waters, bondlen, lat.cell, grid, distance=False)]
            if logger.level <= logging.DEBUG:
                pairs2 = [v for v in pl.pairlist_crude(lat.waters, bondlen, lat.cell, distance=False)]
                logger.debug("pairs: {0}".format(len(pairs)))
                logger.debug("pairs2: {0}".format(len(pairs2)))
                for pair in pairs:
                    i,j = pair
                    assert (i,j) in pairs2 or (j,i) in pairs2
                for pair in pairs2:
                    i,j = pair
                    assert (i,j) in pairs or (j,i) in pairs
            logger.info("  End estimating the bonds.")

        shuffled_pairs = []
        for pair in pairs:
            i,j = pair
            if random.randint(0,1) == 0:
                i,j = j,i
            shuffled_pairs.append((i,j))
            
        graph = dg.IceGraph()
        graph.register_pairs(shuffled_pairs)
        logger.info("End placing the bonds.")

    
    #replicate water molecules to make a repeated cell
    reppositions = replicate(lat.waters, rep)


    #scale the cell
    for d in range(3):
        lat.cell[:,d] *= rep[d]
        
    result = {"positions"   : reppositions,
              "cell"        : lat.cell,
              "celltype"    : lat.celltype,
              "bondlen"     : bondlen}
    if stage == 0:
        return result
        #return reppositions, None, None, lat.cell, lat.celltype, bondlen

    #replicate the graph
    #This also shuffles the bond directions
    graph = replicate_graph(graph, lat.waters, rep)

    result["graph"] = graph
    if stage == 1:
        return result

    #Test
    undir = graph.to_undirected()
    for node in range(undir.number_of_nodes()):
        if node not in undir:
            logger.debug("z=0 at {0}".format(node))
        else:
            z = len(undir.neighbors(node))
            if  z!= 4:
                logger.debug("z={0} at {1}".format(z,node))

    if graph.number_of_edges() != len(reppositions)*2:
        logger.info("Inconsistent number of HBs {0} for number of molecules {1}.".format(graph.number_of_edges(),len(reppositions)))
        return result

    #make them obey the ice rule
    logger.info("Start making the bonds obey the ice rules.")
    graph.purge_ice_defects()
    logger.info("End making the bonds obey the ice rules.")
    if stage == 2:
        result["graph"] = graph
        return result

    #Rearrange HBs to purge the total dipole moment.
    logger.info("Start depolarization.")
    double_net_test = True
    try:
        if lat.double_network:
            double_net_test = (rep[0] % 2 == 0) and (rep[1] % 2 == 0) and (rep[2] % 2 == 0)
    except:
        pass #just ignore.
    if not double_net_test:
        logger.error("In making the ice structure having the double network (e.g. ices 6 and 7), all the repetition numbers (--rep) must be even.")
        sys.exit(1)
    spacegraph = dg.SpaceIceGraph(graph,coord=reppositions)
    draw = dg.YaplotDraw(reppositions, lat.cell, data=spacegraph)
    yapresult  = dg.depolarize(spacegraph, lat.cell, draw=draw)
    logger.info("End depolarization.")
    #determine the orientations of the water molecules based on edge directions.
    rotmatrices = orientations(reppositions, spacegraph, lat.cell)
    result["rotmatrices"] = rotmatrices
    result["graph"]       = spacegraph
    result["yaplot"]      = yapresult
    if stage == 3:
        return result

        




def generate_cages(lattice_type, rep):
    logger = logging.getLogger()
    logger.info("Ice type: {0}".format(lattice_type))
    lat = safe_import("lattice", lattice_type)
    try:
        cagepos, cagetype = parse_cages(lat.cages)
        return replicate(cagepos,      rep), cagetype
    except AttributeError:
        return None, None

