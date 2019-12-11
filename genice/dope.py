import networkx as nx
import sys
import itertools as it
from logging import getLogger, DEBUG, basicConfig
import heapq

def load_NGPH(file):
    g = nx.DiGraph()
    while True:
        line=file.readline()
        if len(line) == 0:
            break
        if len(line)>5 and line[:5] == "@NGPH":
            N = int(file.readline())
            while True:
                x,y = [int(x) for x in file.readline().split()]
                if x < 0:
                    break
                g.add_edge(x,y)
    return g




def shortest_path(G, start, end, occupied_edges, forbidden_vertices):
    """
    Find a shortest path from the start to the end.

    Returns:
    a list of vertices from the start to the end.
    """
    logger = getLogger()
    q = [(0, [start,])]  # Heap of (cost, path)
    visited = set()
    while len(q):
        # logger.debug(q)
        (cost, path) = heapq.heappop(q)
        v0 = path[-1]
        visited.add(v0)
        for v1 in G[v0]:
            if (v0,v1) in occupied_edges:
                continue
            if v1 == end:
                return path+[v1]
            if v1 in forbidden_vertices:
                continue
            if v1 not in visited:
                heapq.heappush(q, (cost + 1, path+[v1]))
    return None
    





def bipartile_self_avoiding_shortest_path_group(g, anions, cations):
    """
    Find a simplest way to embed pairs of ions in the ice lattice.

    Return a set of paths that are
    1. connecting one of anions to one of cations twice.
    2. self-avoiding, i.e. any two paths does not share an edge.
    3. total path length is shortest.
    """
    
    cheapest=999999
    best=None
    forbidden_vertices = frozenset(anions+cations)

    processed = set()
    for cat1 in it.permutations(cations):
        for cat2 in it.permutations(cations):
            occupied_edges = set()
            paths = []
            cost = 0
            ac = tuple([(a,c) for a,c in zip(anions+anions, cat1+cat2)])
            if ac not in processed:
                processed.add(ac)
                #print(ac)
                for a,c in ac:
                    path = shortest_path(g, a, c, occupied_edges, forbidden_vertices)
                    #print(path)
                    #print(a,c)
                    if path is None:
                        # This may happen when ions are too close each other. Some ion may block the path.
                        cost = 999999
                    else:
                        paths.append(path)
                        cost += len(path)-1
                        for k in range(len(path)-1):
                            occupied_edges.add((path[k], path[k+1]))
            if cost < cheapest:
                best = paths
                cheapest = cost
                print(cost)
    return best

if __name__ == "__main__":
    basicConfig(level=DEBUG, format='%(asctime)s- %(name)s - %(levelname)s - %(message)s')
    logger = getLogger(__name__)
    logger.setLevel(DEBUG)
    g = load_NGPH(sys.stdin)
    anions = [0,5,10]
    cations = [4,20,25]
    print(bipartile_self_avoiding_shortest_twin_paths(g, anions, cations))
