import itertools as it

def triangles(neighbors):
    done = set()
    for i in neighbors:
        for j in neighbors[i]:
            for k in neighbors[j]:
                if k in neighbors[i]:
                    s = frozenset((i,j,k))
                    if s not in done:
                        yield s
                    done.add(s)
