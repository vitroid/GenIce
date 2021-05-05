# coding: utf-8

import genice2.formats
from genice2.decorators import timeit, banner
from math import atan2, pi, degrees
from openpyscad import *
from logging import getLogger
import numpy as np
desc = {"ref": {},
        "brief": "OpenSCAD.",
        "usage": """
Usage: genice2 icename -f openscad[options]

Options:
    scale=50
    rnode=0.07
    rbond=0.06
    fn=20
"""
        }


# from cycles.cycles import cycles_iter


# primitives
def rhomb(cell):
    origin = np.zeros(3)
    points = [x + y + z for x in (origin, cell[0])
              for y in (origin, cell[1]) for z in (origin, cell[2])]
    faces = [[0, 1, 3, 2], [0, 4, 5, 1], [0, 2, 6, 4],
             [5, 4, 6, 7], [6, 2, 3, 7], [3, 1, 5, 7]]
    return Polyhedron(points=[list(point) for point in points],
                      faces=[face for face in faces])

# it becomes much easier if the openpyscad supports multmatrix.
# def rhomb(cell):
#    affine = np.zeros([3,4])
#    affine[:3,:3] = cell.T
#    return Cube().multmatrix(affine.tolist())

# copied from rings.py


def bond(p1, p2, r):
    d = p2 - p1
    H = np.linalg.norm(d)
    R = np.linalg.norm(d[:2])
    theta = pi / 2 - atan2(d[2], R)
    phi = atan2(d[1], d[0])
    return Cylinder(r=r, h=H).rotate(
        [0, degrees(theta), degrees(phi)]).translate(list(p1))


def test():
    print(Sphere(r=5).translate([1, 2, 3]).dumps())
    print((Sphere(r=2) + Sphere(r=3)).dumps())


class Format(genice2.formats.Format):
    """
The positions of the water molecules are output in OpenSCAD format.

Options:
  scale=50    Scaling factor.
  rnode=0.07  Radius of the nodes.
  rbond=0.06  Radius of the bonds.
  fn=20       Number of divisions for a curved surface.
  # cycle=False Draw a polygon at each cycle in the HB network.
    """

    def __init__(self, **kwargs):
        self.options = dict(scale=50, rnode=0.07,
                            rbond=0.06, fn=20, cycle=False)
        unknown = dict()
        for k, v in kwargs.items():
            if k in ("scale", "rnode", "rbond", 'fn'):
                self.options[k] = float(v)
            # elif k == "cycle":
            #     self.options[k] = bool(v)
            else:
                unknown[k] = v
        super().__init__(**unknown)

    def hooks(self):
        return {0: self.Hook0, 2: self.Hook2}

    @timeit
    @banner
    def Hook0(self, ice):
        "Preprocess for OpenSCAD."
        for d in range(3):
            ice.rep[d] += 2  # Extend the size,then cut off later.

    @timeit
    @banner
    def Hook2(self, ice):
        "Draw network with OpenSCAD."
        logger = getLogger()
        scale = self.options["scale"]
        rnode = self.options["rnode"]
        rbond = self.options["rbond"]
        fn = self.options["fn"]
        # cycle = self.options["cycle"]
        cellmat = ice.repcell.mat
        rep = np.array(ice.rep)
        trimbox = ice.cell.mat * np.array([(rep[i] - 2) for i in range(3)])
        trimoffset = ice.cell.mat[0] + ice.cell.mat[1] + ice.cell.mat[2]
        # logger.info(ice.repcell.mat)
        # logger.info(ice.cell.mat)

        margin = 0.2  # expansion relative to the cell size
        lower = (1.0 - margin) / rep
        upper = (rep - 1.0 + margin) / rep

        bonds = []
        if rbond > 0.0:
            for i, j in ice.graph.edges(data=False):
                s1 = ice.reppositions[i]
                s2 = ice.reppositions[j]
                d = s2 - s1
                d -= np.floor(d + 0.5)
                logger.debug("Len {0}-{1}={2}".format(i, j, np.linalg.norm(d)))
                s2 = s1 + d
                if ((lower[0] < s1[0] < upper[0] and
                     lower[1] < s1[1] < upper[1] and
                     lower[2] < s1[2] < upper[2]) or
                    (lower[0] < s2[0] < upper[0] and
                     lower[1] < s2[1] < upper[1] and
                     lower[2] < s2[2] < upper[2])):
                    bonds.append((s1 @ cellmat, s2 @ cellmat))

        nodes = []
        if rnode > 0.0:
            for s1 in ice.reppositions:
                if lower[0] < s1[0] < upper[0] and lower[1] < s1[1] < upper[1] and lower[2] < s1[2] < upper[2]:
                    nodes.append(s1 @ cellmat)

        objs = Union()
        for node in nodes:
            objs.append(Sphere(r="Rnode").translate(list(node)))
        for s1, s2 in bonds:
            objs.append(bond(s1, s2, r="Rbond"))
        i = Intersection()
        i.append(objs)
        i.append(rhomb(trimbox).translate(list(trimoffset)))
        s = f"$fn={fn};Rnode={rnode};Rbond={rbond};\n"
        # 文法的に間違ってないようだが、表示されない。
        # Renderボタンを押すと時間がかかるので、
        # ちゃんとRenderすれば表示される、はず。
        self.output = s + \
            i.translate(list(-trimoffset)).scale([scale, scale, scale]).dumps()


if __name__ == "__main__":
    test()
