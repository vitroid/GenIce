"""
Helper functions to read various value types.
"""

import numpy as np

# workaround for python 3.10 and newer
try:
    from collections.abc import Iterable
except ImportError:
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


def parse_guest(guests, arg):
    """
    Parameters
        guests: a dictionary; key is the cage type, value is the guest fractions
        arg: Argumets of --guest option.

    Returns:
        guests: the updated dictionary.
    """
    cagetype, spec = arg.split("=")

    # spec contains a formula consisting of "+" and "*"
    contents = spec.split("+")

    assert cagetype not in guests, "Cage type already specified."

    for content in contents:
        if "*" in content:
            molec, frac = content.split("*")
            frac = float(frac)
            guests[cagetype][molec] = frac
        else:
            molec = content
            frac = 1.0
            guests[cagetype][molec] = frac

    return guests


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


def plugin_option_parser(s):
    """
    Separate the plugin name and options
    """

    left = s.find("[")
    right = s.find("]")
    if 0 < left < len(s) and 0 < right < len(s) and left < right:
        args = s[left + 1 : right]
        name = s[:left]
    else:
        return s, {}

    kwargs = dict()
    for elem in args.split(":"):
        if "=" in elem:
            k, v = elem.split("=", 1)
            kwargs[k] = v
        else:
            kwargs[elem] = True
    return name, kwargs
