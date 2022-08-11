

import string
from logging import getLogger

import networkx as nx
import numpy as np
from cycless.cycles import centerOfMass, cycles_iter
from cycless.polyhed import cage_to_graph, polyhedra_iter
# for cage assessment
from graphstat import GraphStat


def assess_cages(graph, node_pos):
    """Assess cages from  the graph topology.

    Args:
        graph (graph-like): HB network
        nodepos (np.Array): Positions of the nodes
    """
    logger = getLogger()

    # Prepare the list of rings
    # taking the positions in PBC into account.
    ringlist = [[int(x) for x in ring] for ring in cycles_iter(nx.Graph(graph), 8, pos=node_pos)]

    # Positions of the centers of the rings.
    ringpos = [centerOfMass(ringnodes, node_pos) for ringnodes in ringlist]

    maxcagesize = 22
    cagepos = []
    cagetypes = []
    # data storage of the found cages
    db = GraphStat()
    labels = set()
    g_id2label = dict()

    # Detect cages and classify
    for cage in polyhedra_iter(ringlist, maxcagesize):
        cagepos.append(centerOfMass(list(cage), ringpos))
        g = cage_to_graph(cage, ringlist)
        cagesize = len(cage)
        g_id = db.query_id(g)
        if g_id < 0:
            g_id = db.register()
            enum = 0
            label = f"A{cagesize}"
            while label in labels:
                char = string.ascii_lowercase[enum]
                label = f"A{cagesize}{char}"
                enum += 1
            g_id2label[g_id] = label
            labels.add(label)
            ringcount = [0 for i in range(9)]
            for ring in cage:
                ringcount[len(ringlist[ring])] += 1
            index = []
            for i in range(9):
                if ringcount[i] > 0:
                    index.append(f"{i}^{ringcount[i]}")
            index = " ".join(index)
            logger.info(f"    Cage type: {label} ({index})")

        else:
            label = g_id2label[g_id]
        cagetypes.append(label)
    if len(cagepos) == 0:
        logger.info("    No cages detected.")
    return np.array(cagepos), cagetypes
