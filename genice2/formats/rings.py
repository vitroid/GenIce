# coding: utf-8
"""
Show rings in Yaplot format, defined in https://github.com/vitroid/Yaplot

Usage:
    genice2 III -f rings > 3.rings.yap
    genice2 III -f rings[max=5] > output.yap    # up to 5-membered rings.
    genice2 III -f rings[openscad:max=5] > output.scad    # up to 5-membered rings.
"""

from genice2.decorators import timeit, banner
from cycless.cycles import cycles_iter
import genice2.formats
from genice2 import rigid
import colorsys
from openpyscad import *
from logging import getLogger
import yaplotlib as yp
import networkx as nx
import numpy as np
from math import atan2, pi, degrees
from collections import defaultdict
import sys
desc = {
    "ref": {
        "MBO2007": "Matsumoto, M., Baba, A. & Ohmine, I. Topological building blocks of hydrogen bond network in water. J. Chem. Phys. 127, 134504 (2007).",
        "CountRings": "https://github.com/vitroid/CountRings",
        "Yaplot": "https://github.com/vitroid/Yaplot"},
    "brief": "Show rings in Yaplot.",
    "usage": __doc__}


# Basic primitives
# Store everything in an array
class Prims:
    def __init__(self):
        self.prims = []
        self.palette = {}

    def SetPalette(self, num, color):
        self.palette[num] = color

    def Cylinder(self, p1, p2, r, **kwargs):
        self.prims.append(["cyl", (p1, p2, r), kwargs])

    def Sphere(self, p, r, **kwargs):
        self.prims.append(["sph", (p, r), kwargs])

    def Polygon(self, points, **kwargs):
        self.prims.append(["ply", (points,), kwargs])

    def Render(self):
        for prim in self.prims:
            print(prim)


class YaPrims(Prims):
    @timeit
    @banner
    def Render(self, file=None):
        "Rendering with Yaplot."
        logger = getLogger()
        lastLayer = -1
        lastColor = -1
        s = ""
        for typ, args, kwargs in self.prims:
            if typ == "cyl":
                p1, p2, r = args
                s += yp.Size(r)
                s += yp.Layer(kwargs["layer"])
                s += yp.Color(kwargs["color"])
                s += yp.Line(p1, p2)
            elif typ == "ply":
                points = args[0]
                s += yp.Layer(kwargs["layer"])
                s += yp.Color(kwargs["color"])
                s += yp.Polygon(points)
            else:
                logger.warn("Unknown primitive: {0}".format(typ))
        s += yp.NewPage()
        if file is None:
            return s
        file.write(s)


def bond(p1, p2, r):
    d = p2 - p1
    H = np.linalg.norm(d)
    R = np.linalg.norm(d[:2])
    theta = pi / 2 - atan2(d[2], R)
    phi = atan2(d[1], d[0])
    return Cylinder(r=r, h=H).rotate(
        [0, degrees(theta), degrees(phi)]).translate(list(p1))


class ScadPrims(Prims):
    @timeit
    @banner
    def Render(self, file=None):
        "Rendering with OpenSCAD."
        logger = getLogger()
        layers = defaultdict(Union)
        for typ, args, kwargs in self.prims:
            if typ == "cyl":
                p1, p2, r = args
                layer = kwargs["layer"]
                #pal   = "palette{0}".format(kwargs["color"])
                layers[layer].append(bond(p1, p2, r).color(
                    self.palette[kwargs["color"]]))
            elif typ == "ply":
                p = args[0]
                r = 0.02
                layer = kwargs["layer"]
                #pal   = "palette{0}".format(kwargs["color"])
                for i in range(len(p)):
                    layers[layer].append(Sphere(r).translate(
                        list(p[i])).color(self.palette[kwargs["color"]]))
                    layers[layer].append(
                        bond(p[i - 1], p[i], r).color(self.palette[kwargs["color"]]))
            else:
                logger.warn("Unknown primitive: {0}".format(typ))
        if file is None:
            output = "$fn=20;"
            # for p, color in self.palette.items():
            #    output += "palette{0}={1};\n".format(p,list(color))
            for layer in layers:
                output += "module layer{0}()".format(layer) + \
                    "{" + layers[layer].dumps() + "}\n"
            for layer in layers:
                output += "layer{0}();\n".format(layer)
            return output
        assert False


def face(p, center, rpos):
    pos = rpos * 0.8 + center
    n = rpos.shape[0]
    p.Polygon(pos, color=n, layer=n)


class Format(genice2.formats.Format):
    largestring = 8

    def __init__(self, **kwargs):
        logger = getLogger()
        self.format = "yaplot"
        unknown = dict()
        for k, v in kwargs.items():
            if k == "openscad":
                self.format = k
            elif k == "max":
                self.largestring = int(v)
            else:
                logger.warn("Unknown keyword: {0}".format(k))
        logger.info("  Largest ring: {0}.".format(self.largestring))
        if self.format == "yaplot":
            self.p = YaPrims()
        elif self.format == "openscad":
            self.p = ScadPrims()
        phi = (1 + 5**0.5) / 2  # golden ratio
        self.p.SetPalette(0, [0, 0, 0])
        for i in range(3, self.largestring + 1):
            hue = (i / phi) % 1
            r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 1)
            self.p.SetPalette(i, [r, g, b])
        super().__init__()

    def hooks(self):
        return {2: self.Hook2}

    @timeit
    @banner
    def Hook2(self, ice):
        "Show rings."
        logger = getLogger()
        # copied from svg_poly
        graph = nx.Graph(ice.graph)  # undirected
        cellmat = ice.repcell.mat
        for i, j in graph.edges():
            pi, pj = ice.reppositions[i], ice.reppositions[j]
            d = pj - pi
            d -= np.floor(d + 0.5)
            self.p.Cylinder(pi @ cellmat, (pi + d) @ cellmat, 0.01,  # 0.2 AA
                            color=0,
                            layer=2)
        for ring in cycles_iter(graph, self.largestring, pos=ice.reppositions):
            deltas = np.zeros((len(ring), 3))
            for k, i in enumerate(ring):
                d = ice.reppositions[i] - ice.reppositions[ring[0]]
                d -= np.floor(d + 0.5)
                deltas[k] = d
            comofs = np.sum(deltas, axis=0) / len(ring)
            deltas -= comofs
            com = ice.reppositions[ring[0]] + comofs
            com -= np.floor(com)
            # rel to abs
            com = com @ cellmat
            deltas = deltas @ cellmat
            face(self.p, com, deltas)
        self.output = self.p.Render()
