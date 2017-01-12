#!/usr/bin/env python3
# -*- python -*-

import sys
import numpy     as np
import pairlist  as pl
import digraph   as dg
import random
import math
import logging

def usage(parser):
    parser.print_help()
    sys.exit(1)



def orientations(coord, graph, cell):
    rotmatrices = []
    for node in range(graph.number_of_nodes()):
        nei = graph.neighbors(node)
        oh1 = coord[nei[0]] - coord[node]
        oh1 -= np.floor( oh1 + 0.5 )
        oh1 = np.dot(oh1,cell)                #abs coord
        oh2 = coord[nei[1]] - coord[node]
        oh2 -= np.floor( oh2 + 0.5 )
        oh2 = np.dot(oh2,cell)                #abs coord
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






def generate_ice(lattice_type, density=-1, seed=1000, rep=(1,1,1), noGraph=False, yaplot=False):
    logger = logging.getLogger()
    lat     = __import__(lattice_type)
    lat.waters = np.fromstring(lat.waters, sep=" ")
    lat.waters = lat.waters.reshape((lat.waters.size//3,3))
    #prepare cell transformation matrix
    if lat.celltype == "rect":
        lat.cell = np.fromstring(lat.cell, sep=" ")
        lat.cell = np.diag(lat.cell)
    elif lat.celltype == "monoclinic":
        lat.cell = np.fromstring(lat.cell, sep=" ")
        beta = lat.cell[3] * math.pi / 180.
        lat.cell = np.array(((lat.cell[0]*1.0, lat.cell[1]*0.0, lat.cell[2]*math.cos(beta)),
                             (lat.cell[0]*0.0, lat.cell[1]*1.0, lat.cell[2]*0.0),
                             (lat.cell[0]*0.0, lat.cell[1]*0.0, lat.cell[2]*math.sin(beta))))
        lat.cell = lat.cell.transpose()   #all the vector calculations are done in transposed manner.
    else:
        logger.error("unknown cell type: {0}".format(lat.celltype))
        sys.exit(1)

    #express molecular positions in the coordinate relative to the cell
    if lat.coord == "absolute":
        lat.waters = np.dot(lat.waters,np.linalg.inv(lat.cell),)
        lat.waters = np.array(lat.waters)
    random.seed(seed)

    #Prearranged network topology information (if available)
    pairs = None
    bondlen = None
    try:
        lines = lat.pairs.split("\n")
        pairs = set()
        for line in lines:
            columns = line.split()
            if len(columns) == 2:
                i,j = [int(x) for x in columns]
                if j<i:
                    i,j = j,i
                pairs.add((i,j))
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


    if not noGraph:
        if pairs is None:
            #make bonded pairs according to the pair distance.
            #make before replicating them.
            grid = pl.determine_grid(lat.cell, bondlen)
            pairs = pl.pairlist_fine(lat.waters, bondlen, lat.cell, grid)

        graph = dg.IceGraph()
        graph.register_pairs(pairs)

    
    #replicate water molecules to make a repeated cell
    reppositions = replicate(lat.waters, rep)


    #scale the cell
    for d in range(3):
        lat.cell[:,d] *= rep[d]
        
    if noGraph:
        return reppositions, None, None, lat.cell, lat.celltype, bondlen

    #replicate the graph
    #This also shuffles the bond directions
    graph = replicate_graph(graph, lat.waters, rep)


    #Test
    if graph.number_of_edges() != len(reppositions)*2:
        logger.error("Inconsistent number of HBs.[2] {0} {1}".format(graph.number_of_edges(),len(reppositions)*2))
        sys.exit(1)


    #make them obey the ice rule
    logger.info("Start making the bonds obey the ice rules.")
    graph.purge_ice_defects()
    logger.info("End making the bonds obey the ice rules.")


    #Rearrange HBs to purge the total dipole moment.
    logger.info("Start depolarization.")
    spacegraph = dg.SpaceIceGraph(graph,coord=reppositions)
    #spacegraph.depolarize()
    draw = None
    if yaplot:
        draw = dg.YaplotDraw(reppositions, lat.cell, data=spacegraph)
    yapresult  = dg.depolarize(spacegraph, lat.cell, draw=draw)
    #yapresult  = dg.depolarize(spacegraph,draw=draw)
    logger.info("End depolarization.")


    #determine the orientations of the water molecules based on edge directions.
    rotmatrices = orientations(reppositions, spacegraph, lat.cell)

    result = {"positions"   : reppositions,
              "rotmatrices" : rotmatrices,
              "graph"       : spacegraph,
              "cell"        : lat.cell,
              "celltype"    : lat.celltype,
              "bondlen"     : bondlen,
              "yaplot"      : yapresult}
    return result



def generate_cages(lattice_type, rep):
    lat     = __import__(lattice_type)
    try:
        cagetype = []
        cagepos = []
        for line in lat.cages.split("\n"):
            cols = line.split()
            if len(cols)>0:
                cagetype.append(cols[0])
                cagepos.append([float(x) for x in cols[1:4]])
        cagepos = np.array(cagepos)
        return replicate(cagepos,      rep), cagetype
    except AttributeError:
        return None, None

