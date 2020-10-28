#!/usr/bin/env python3
#
# Simplified from github/Vitrite

import sys
import networkx as nx
from logging import getLogger
from collections import defaultdict
import numpy as np

def cage_to_graph(cage, ringlist):
    g = nx.Graph()
    for ring in cage:
        nodes = ringlist[ring]
        for i in range(len(nodes)):
            g.add_edge(nodes[i-1], nodes[i])
    return g


def centerOfMass(members, rpos):
    logger = getLogger()
    dsum = np.zeros(3)
    for member in members:
        d = rpos[member] - rpos[members[0]]
        d -= np.floor(d+0.5)
        dsum += d
    com = rpos[members[0]] + dsum / len(members)
    com -= np.floor(com)
    return com



#reorder the ring noders so as to start from the first node
def reorder(ring,first,second):
    s = ring.index(first)
    if ring[s-1] == second:
        r = [ring[i] for i in range(s,s-len(ring),-1)]
    else:
        r = [ring[i] for i in range(s-len(ring),s)]
    return r


#get two lists of nodes (rings) and make a large ring.
def MergeRings(ring1,ring2,first,second):
    logger = getLogger()
    r1 = reorder(ring1,first,second)
    r2 = reorder(ring2,first,second)
    logger.debug("#{0}+{1}".format(r1,r2))
    #zipper motion
    head=0
    while r1[head-1] == r2[head-1]:
        head -= 1
        if head == -len(r1):
            return []
    tail=1
    while r1[tail+1-len(r1)] == r2[tail+1-len(r2)]:
        tail += 1
    #unshared nodes of the rings
    rest1 = set(r1) - set([r1[i] for i in range(head,tail+1)])
    rest2 = set(r2) - set([r2[i] for i in range(head,tail+1)])
    #if the remaining parts of the ring have common nodes,
    if len(rest1 & rest2) != 0:
        #not a simple ring
        return None
    ring = [r1[i] for i in range(tail-len(r1),head)]
    ring += [r2[i] for i in range(head,tail-len(r2),-1)]
    logger.debug("#{0} {1} {2} {3} {4}".format(head,tail,ring,rest1,rest2,rest1 & rest2))
    return ring


def Triplets(nodes):
    tri = []
    for i in range(len(nodes)):
        tri.append((nodes[i-2],nodes[i-1],nodes[i]))
    return tri


def Edges(nodes):
    ed = []
    for i in range(len(nodes)):
        ed.append((nodes[i-1],nodes[i]))
    return ed


#Look up all the polyhedral fragments (vitrites) in the given set of rings.
def Polyhed(_rings, maxfragsize=20):
    #Local functions

    def RegisterTriplets(nodes,ringid):
        for triplet in Triplets(nodes):
            _RingsAtATriplet[triplet].append(ringid)
            tr = tuple(reversed(triplet))
            _RingsAtATriplet[tr].append(ringid)

    def RegisterEdges(nodes,ringid):
        for edge in Edges(nodes):
            _RingsAtAnEdge[edge].append(ringid)
            ed = tuple(reversed(edge))
            _RingsAtAnEdge[ed].append(ringid)

    def IsDivided(fragment):
        nodes = set()
        for ring in fragment:
            nodes |= set(_rings[ring])
        G2 = _G.copy()
        ebunch = []
        for i in nodes:
            for j in _G.neighbors(i):
                ebunch.append((i,j))
        #G2.remove_edges_from(ebunch)
        G2.remove_nodes_from(nodes)
        logger.debug("NCOMPO: {0} {1}".format(nx.number_connected_components(G2),_ncompo))
        return nx.number_connected_components(G2) > _ncompo

    #Return True if the given fragment contains rings that are not the member of the fragment.
    def ContainsExtraRing(fragment):
        tris = set()
        allnodes = set()
        #A fragment is a set of ring IDs.
        for ringid in fragment:
            nodes = _rings[ringid]
            allnodes |= set(nodes)
            tris |= set(Triplets(nodes))
        for tri in tris:
            for ringid in _RingsAtATriplet[tri]:
                if ringid not in fragment:
                    #if all the nodes of a ring is included in the fragment,
                    nodes = _rings[ringid]
                    if len(set(nodes) - allnodes) == 0:
                        return True
                    #logger.debug(fragment,ringid)
        return False
    #Grow up a set of rings by adding new ring on the perimeter.
    #Return the list of polyhedron.
    #origin is the ring ID of the first face in the polyhedron
    def Progress(origin, peri, fragment, numRingsOnTheNode):
        #Here we define a "face" as a ring belonging to the (growing) polyhedron.
        logger.debug("#{0} {1}".format(peri,fragment))
        if len(fragment) > maxfragsize:
            #logger.debug("#LIMITTER")
            return False
        #if the perimeter closes,
        if len(peri) == 0:
            #If the polyhedron has internal vertices that are not a part of the polyhedron (i.e. if the polyhedron does not divide the total network into to components)
            if not IsDivided(fragment):
                #If the fragment does not contain any extra ring whose all vertices belong to the fragment but the ring is not a face,
                if not ContainsExtraRing(fragment):
                    #Add the fragment to the list.
                    #A fragment is a set of ring IDs of the faces
                    _vitrites.add(frozenset(fragment))
                else:
                    logger.debug("It contains extra rings(s).")
            else:
                logger.debug("It has internal vertices.")
            #Search finished.
            return True
        #If the perimeter is still open,
        for i in range(len(peri)):
            #If any vertex on the perimeter is shared by more than two faces,
            if numRingsOnTheNode[peri[i]] > 2:
                logger.debug("#Failed(2)")
                return False
        for i in range(len(peri)):
            #Look up the node on the perimeter which is shared by two faces.
            if numRingsOnTheNode[peri[i]] == 2:
                #Reset the frag
                trynext = False
                #Three successive nodes around the node i
                center = peri[i]
                left   = peri[i-1]
                right  = peri[i+1-len(peri)] #Avoid to refer the out-of-list element
                #Reset the frag
                success = False
                logger.debug("Next triplet:{0} {1} {2}".format(left,center,right))
                if (left,center,right) in _RingsAtATriplet:
                    logger.debug("Here rings are:{0}".format(_RingsAtATriplet[(left,center,right)]))
                    for ringid in _RingsAtATriplet[(left,center,right)]:
                        logger.debug("#Next:{0}".format(ringid))
                        #if the ring is new and its ID is larger than the origin,
                        if origin < ringid and not ringid in fragment:
                            nodes = _rings[ringid]
                            #Add the ring as a face and extend the perimeter
                            newperi = MergeRings(peri, nodes, center,right)
                            logger.debug("#Result:{0}".format(newperi))
                            # result is not a simple ring
                            if newperi == None:
                                trynext = True
                                logger.debug("#Try next!")
                            else:
                                for node in nodes:
                                    numRingsOnTheNode[node] += 1
                                    mult = [numRingsOnTheNode[i] for i in newperi]
                                    logger.debug("#{0} {1} {2} {3} {4}".format(peri, nodes, edge, newperi,mult))
                                result = Progress(origin, newperi, fragment | set([ringid,]), numRingsOnTheNode)
                                for node in nodes:
                                    numRingsOnTheNode[node] -= 1
                                #if result == True:
                                #    return True
                #it might be too aggressive
                if not trynext:
                    break
        logger.debug("#Failed to expand perimeter {0} {1}".format(peri,fragment))
        return False

    logger = getLogger()

    _RingsAtATriplet = defaultdict(list)
    _RingsAtAnEdge = defaultdict(list)

    for ringid, ring in enumerate(_rings):
        RegisterTriplets(ring,ringid)
        RegisterEdges(ring,ringid)
    #For counting the number of components separated by a polyhedral fragment
    _G = nx.Graph()

    for ring in _rings:
        nx.add_cycle(_G,ring)
    _ncompo = nx.number_connected_components(_G)

    _vitrites = set()
    #The first ring
    for ringid in range(len(_rings)):
        peri = _rings[ringid]
        fragment = set([ringid])
        edge = tuple(peri[0:2])
        numRingsOnTheNode = defaultdict(int)
        #increment the number-of-rings-at-a-node counter
        #for each node on the first ring.
        for node in peri:
            numRingsOnTheNode[node] = 1
            logger.debug("#Candid: {0} {1}".format(ringid,_RingsAtAnEdge[edge]))
        #The second ring, which is adjacent to the first one.
        for ringid2 in _RingsAtAnEdge[edge]:
            #The second one must have larger ring ID than the first one.
            if ringid < ringid2:
                nodes = _rings[ringid2]
                #Make the perimeter of two rings.
                newperi = MergeRings(peri, nodes, edge[0],edge[1])
                if newperi != None:
                    #increment the number-of-rings-at-a-node counter
                    #for each node on the second ring.
                    for node in nodes:
                        numRingsOnTheNode[node] += 1
                        mult = [numRingsOnTheNode[i] for i in newperi]
                        logger.debug("{0} {1} {2} {3} {4}".format(peri, nodes, edge, newperi,mult))
                    #Expand the perimeter by adding new faces to the polyhedron.
                    Progress(ringid, newperi, set([ringid,ringid2]), numRingsOnTheNode)
                    #decrement the number-of-rings-at-a-node counter
                    #for each node on the second ring.
                    for node in nodes:
                        numRingsOnTheNode[node] -= 1
    return _vitrites
