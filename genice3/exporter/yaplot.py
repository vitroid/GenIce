# coding: utf-8

import itertools as it
import sys
from collections import defaultdict
from logging import getLogger
from io import TextIOWrapper

import numpy as np
import yaplotlib as yp

from genice3.genice import GenIce3
from genice3.exporter import (
    parse_guest_option,
    parse_spot_guest_option,
    parse_water_model_option,
)
from genice3.util import serialize

desc = {
    "ref": {"Codes": "https://github.com/vitroid/Yaplot"},
    "brief": "Yaplot.",
    "usage": """
Usage: genice2 icename -f yaplot[options]

options:
    H=x   Set the radius of H to be x.
""",
}


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **options):
    "Draw the cell in Yaplot format."
    logger = getLogger()

    size_O = 0.03
    size_H = float(options.get("H", 0)) * size_O

    s = yp.Layer(2)
    x, y, z = genice.cell
    for p, q, r in ((x, y, z), (y, z, x), (z, x, y)):
        for a in (np.zeros(3), p, q, p + q):
            s += yp.Line(a, a + r)

    if size_H == 0:
        # 水素のサイズ指定がない場合
        # prepare the reverse dict
        waters = defaultdict(dict)
        pos = genice.lattice_sites
        s += yp.Layer(4)
        s += yp.Color(3)
        s += yp.Size(0.03)
        for p in pos:
            s += yp.Circle(p @ genice.cell)
        s += yp.Layer(5)
        s += yp.Color(4)
        # s += yp.Size(0.03)
        for i, j in genice.graph.edges(data=False):
            # if i in waters and j in waters:  # edge may connect to the dopant
            O1, O2 = pos[i], pos[j]
            d = O2 - O1
            d -= np.floor(d + 0.5)
            O2 = O1 + d
            s += yp.Line(O1 @ genice.cell, O2 @ genice.cell)

    # 設定可能なオプションはguestとspot_guest。
    guest_info = parse_guest_option(options.get("guest", {}))
    spot_guest_info = parse_spot_guest_option(options.get("spot_guest", {}))
    # waterとwater_modelの両方をサポート（後方互換性のため）
    water_model_name = options.get("water_model") or options.get("water", "4site")
    water_model = parse_water_model_option(water_model_name)
    # water = FourSiteWater()  # dummy

    waters = genice.water_molecules(water_model=water_model)
    guests = genice.guest_molecules(guests=guest_info, spot_guests=spot_guest_info)
    ions = genice.substitutional_ions()

    atoms = serialize(list(waters.values()))

    logger.info("  Total number of atoms: {0}".format(len(atoms)))

    # prepare the reverse dict
    water_sites = defaultdict(dict)
    for water_index, water in waters.items():
        for atom_name, position in zip(water.labels, water.sites):
            if "O" in atom_name:
                water_sites[water_index]["O"] = position
            elif "H" in atom_name:
                if "H0" not in water_sites[water_index]:
                    water_sites[water_index]["H0"] = position
                else:
                    water_sites[water_index]["H1"] = position
    s += yp.Color(3)
    for water_index, water in water_sites.items():
        O = water["O"]
        H0 = water["H0"]
        H1 = water["H1"]
        s += yp.Layer(4)
        s += yp.Color(3)
        s += yp.Size(size_O)
        s += yp.Circle(O)

        s += yp.Line(O, H0)
        s += yp.Line(O, H1)

        s += yp.Size(size_H)
        s += yp.Circle(H0)
        s += yp.Circle(H1)

        s += yp.Color(2)
        s += yp.Layer(1)
        s += yp.Text(O, f"{water_index}")
    s += yp.Layer(3)
    s += yp.Color(4)
    s += yp.ArrowType(1)
    s += yp.Size(size_O)
    for i, j in genice.digraph.edges(data=False):
        if i in waters and j in waters:  # edge may connect to the dopant
            O = water_sites[j]["O"]
            H0 = water_sites[i]["H0"]
            H1 = water_sites[i]["H1"]
            d0 = H0 - O
            d1 = H1 - O
            rr0 = d0 @ d0
            rr1 = d1 @ d1
            if rr0 < rr1 and rr0 < 0.245**2:
                s += yp.Arrow(H0, O)
            if rr1 < rr0 and rr1 < 0.245**2:
                s += yp.Arrow(H1, O)

    gatoms = serialize(guests) + serialize(list(ions.values()))
    palettes = dict()

    s += yp.Layer(4)
    s += yp.ArrowType(1)
    H = []
    O = ""
    for molecule_name, atom_name, position in gatoms:
        if atom_name in palettes:
            pal = palettes[atom_name]
        else:
            pal = 4 + len(palettes)
            palettes[atom_name] = pal
        s += yp.Color(pal)
        if atom_name[0] == "H":
            s += yp.Size(size_H)
        else:
            s += yp.Size(size_O)
        s += yp.Circle(position)
    for a, b in it.combinations(gatoms, 2):
        resname, atomname, position1 = a
        resname, atomname, position2 = b
        d = position1 - position2
        if d @ d < 0.16**2:
            s += yp.Line(position1, position2)
    s += "#" + "\n#Command line: " + " ".join(sys.argv) + "\n"
    s += yp.NewPage()
    file.write(s)
