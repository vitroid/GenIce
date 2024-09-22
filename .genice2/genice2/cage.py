import string
from logging import getLogger

import networkx as nx
import numpy as np
from cycless.cycles import centerOfMass, cycles_iter
from cycless.polyhed import cage_to_graph, polyhedra_iter

# for cage assessment
from graphstat import GraphStat


def assign_unused_label(basename, labels):
    enum = 0
    label = f"A{basename}"
    while label in labels:
        char = string.ascii_lowercase[enum]
        label = f"A{basename}{char}"
        enum += 1
    return label


def make_cage_expression(ring_ids, ringlist):
    ringcount = [0 for i in range(9)]
    for ring in ring_ids:
        ringcount[len(ringlist[ring])] += 1
    index = []
    for i in range(9):
        if ringcount[i] > 0:
            index.append(f"{i}^{ringcount[i]}")
    index = " ".join(index)
    return index


def assess_cages(graph, node_pos):
    """Assess cages from  the graph topology.

    Args:
        graph (graph-like): HB network
        nodepos (np.Array): Positions of the nodes
    """
    logger = getLogger()

    # Prepare the list of rings
    # taking the positions in PBC into account.
    ringlist = [
        [int(x) for x in ring] for ring in cycles_iter(nx.Graph(graph), 8, pos=node_pos)
    ]

    # Positions of the centers of the rings.
    ringpos = [centerOfMass(ringnodes, node_pos) for ringnodes in ringlist]

    MaxCageSize = 22
    cagepos = []
    cagetypes = []
    # data storage of the found cages
    db = GraphStat()
    labels = set()
    g_id2label = dict()

    # Detect cages and classify
    cages = [cage for cage in polyhedra_iter(ringlist, MaxCageSize)]
    cagepos = [centerOfMass(list(cage), ringpos) for cage in cages]
    for cage in cages:
        g = cage_to_graph(cage, ringlist)
        cagesize = len(cage)
        g_id = db.query_id(g)
        # if it is a new cage type
        if g_id < 0:
            # new type!
            # register the last query
            g_id = db.register()

            # prepare a new label
            label = assign_unused_label(cagesize, labels)
            g_id2label[g_id] = label
            labels.add(label)

            # cage expression
            index = make_cage_expression(cage, ringlist)
            logger.info(f"    Cage type: {label} ({index})")
        else:
            label = g_id2label[g_id]
        cagetypes.append(label)
    if len(cagepos) == 0:
        logger.info("    No cages detected.")
    return np.array(cagepos), cagetypes
