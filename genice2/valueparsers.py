import numpy as np
import logging
from collections import Iterable


def flatten(item):
    """Yield items from any nested iterable; see REF."""
    if isinstance(item, str):
        yield item
    elif isinstance(item, Iterable):
        for x in item:
            yield from flatten(x)
    else:
        yield item


def parse_cages(cages):
    cagetype = []
    cagepos = []
    if isinstance(cages, str):
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
    cagepos -= np.floor(cagepos)
    return cagepos, cagetype


def parse_pairs(values):
    if isinstance(values, str):
        lines = values.split("\n")
        pairs = []
        for line in lines:
            columns = line.split()
            if len(columns) == 2:
                i, j = [int(x) for x in columns]
                pairs.append((i, j))
        return pairs
    elif isinstance(values, list):
        pairs = []
        for pair in values:
            pairs.append(tuple(pair[:2]))  # Make it a tuple
        return pairs


def put_in_array(v):
    """
    Obtain a Numpy array from any kind of values.
    """
    if isinstance(v, str):
        return np.fromstring(v, sep=" ")
    elif isinstance(v, list):
        return np.array(v)
    else:
        return v
