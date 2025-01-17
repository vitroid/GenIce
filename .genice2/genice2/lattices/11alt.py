#!/usr/bin/env python
"""
c軸方向に分極した強誘電体氷11の単位胞を、a軸方向に上向き(u)または下向き(d)に並べ、a軸方向に長い単位胞を生成する。
uとdを切りかえるために、間に自動的に遷移層が挿入される。

例えば、genice2 11alt[uuuddd] と指定すると、実際にはuとdの間に分極のない遷移層"_"をはさんだ_uuu_dddという形の結晶が生成する。
"""
import genice2.lattices
import numpy as np
from genice2.cell import cellvectors
from logging import getLogger

desc = {
    "ref": {},
    "usage": "No options available.",
    "brief": "A layered ferroelectric Ice XI.",
    "test": ({"args": "udud"},),
}


class Lattice(genice2.lattices.Lattice):
    def __init__(self, **kwargs):
        logger = getLogger()

        sites = np.array(
            [
                [1 / 2, 0, 1 / 8],
                [1 / 2, 1 / 3, 0],
                [0, 1 / 2, 1 / 8],
                [0, 5 / 6, 0],
                [1 / 2, 0, 1 / 2],
                [1 / 2, 1 / 3, 5 / 8],
                [0, 1 / 2, 1 / 2],
                [0, 5 / 6, 5 / 8],
            ]
        )
        # c軸方向に上向きまたは下向きに分極した単位胞をa軸方向に重ねる。

        _layers = "uu"

        # parse options
        for k, v in kwargs.items():
            if k == "layers":
                _layers = v
            elif v is True:
                # unlabeled option
                _layers = k

        last = _layers[-1]
        layers = ""
        for layer in _layers:
            if last != layer:
                if layer == "u":
                    layers += "/"
                else:
                    layers += "\\"
            layers += layer
            last = layer
        logger.info(f"{layers}; {len(layers)} layers along the 'a' axis.")

        N_layers = len(layers)

        assert N_layers > 1, "More than one layer must be specified."

        edges = []
        waters = []
        last = layers[-1]
        for i, layer in enumerate(layers):
            waters.append(sites + np.array([i, 0, 0]))
            offset = i * 8
            right = (i + 1) * 8
            if i + 1 == N_layers:
                right = 0
            if layer == "u":
                assert last in ("u", "/"), "Incompatible u layer."
                edges.append([offset + 0, offset + 1])
                edges.append([offset + 1, offset + 2])
                edges.append([offset + 2, offset + 3])
                edges.append([offset + 3, offset + 0])
                edges.append([offset + 1, right + 2])
                edges.append([right + 3, offset + 0])
                edges.append([offset + 0, offset + 4])
                edges.append([offset + 2, offset + 6])
                edges.append([offset + 4, offset + 7])
                edges.append([offset + 7, offset + 6])
                edges.append([offset + 6, offset + 5])
                edges.append([offset + 5, offset + 4])
                edges.append([offset + 4, right + 7])
                edges.append([right + 6, offset + 5])
                edges.append([offset + 5, offset + 1])
                edges.append([offset + 7, offset + 3])
            elif layer == "d":
                assert last in ("d", "\\"), "Incompatible d layer."
                edges.append([offset + 0, offset + 3])
                edges.append([offset + 3, offset + 2])
                edges.append([offset + 2, offset + 1])
                edges.append([offset + 1, offset + 0])
                edges.append([right + 2, offset + 1])
                edges.append([offset + 0, right + 3])
                edges.append([offset + 4, offset + 0])
                edges.append([offset + 6, offset + 2])
                edges.append([offset + 7, offset + 4])
                edges.append([offset + 6, offset + 7])
                edges.append([offset + 5, offset + 6])
                edges.append([offset + 4, offset + 5])
                edges.append([right + 7, offset + 4])
                edges.append([offset + 5, right + 6])
                edges.append([offset + 1, offset + 5])
                edges.append([offset + 3, offset + 7])
            elif layer == "/":
                assert last in ("d", "\\"), "Incompatible d2u layer."
                edges.append([offset + 0, offset + 3])
                edges.append([offset + 3, offset + 2])
                edges.append([offset + 2, offset + 1])
                edges.append([offset + 0, offset + 1])
                edges.append([offset + 1, right + 2])
                edges.append([right + 3, offset + 0])
                edges.append([offset + 4, offset + 0])
                edges.append([offset + 6, offset + 2])

                edges.append([offset + 5, offset + 4])
                edges.append([offset + 5, offset + 6])
                edges.append([offset + 6, offset + 7])
                edges.append([offset + 7, offset + 4])
                edges.append([offset + 4, right + 7])
                edges.append([right + 6, offset + 5])
                edges.append([offset + 1, offset + 5])
                edges.append([offset + 3, offset + 7])
            elif layer == "\\":
                assert last in ("u", "/"), "Incompatible u2d layer."
                edges.append([offset + 1, offset + 0])
                edges.append([offset + 1, offset + 2])
                edges.append([offset + 2, offset + 3])
                edges.append([offset + 3, offset + 0])
                edges.append([right + 2, offset + 1])
                edges.append([offset + 0, right + 3])
                edges.append([offset + 0, offset + 4])
                edges.append([offset + 2, offset + 6])

                edges.append([offset + 4, offset + 7])
                edges.append([offset + 7, offset + 6])
                edges.append([offset + 6, offset + 5])
                edges.append([offset + 4, offset + 5])
                edges.append([right + 7, offset + 4])
                edges.append([offset + 5, right + 6])
                edges.append([offset + 5, offset + 1])
                edges.append([offset + 7, offset + 3])
            last = layer

        self.waters = np.vstack(waters) / np.array([N_layers, 1, 1])
        N_waters = N_layers * 8
        logger.debug(self.waters.shape)
        self.fixed = edges

        self.pairs = self.fixed

        a = 4.4923 / 10 * N_layers
        b = 7.7808 / 10
        c = 7.3358 / 10

        self.cell = cellvectors(a, b, c)

        self.density = 18 * N_waters / 6.022e23 / (np.linalg.det(self.cell) * 1e-21)

        self.coord = "relative"
