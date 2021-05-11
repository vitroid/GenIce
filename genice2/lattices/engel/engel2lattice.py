#!/usr/bin/env python3
# coding: utf-8
import os
import json
import numpy as np


def load_castep(file):
    mode = ""
    for line in file:
        if line[0] == "%":
            if line[1:6] == "BLOCK":
                cols = line.split()
                if cols[1][:3] == "pos":
                    mode = "pos"
                    Hs = []
                    Os = []
                elif cols[1][:3] == "lat":
                    mode = "lat"
                    cell = []
            elif line[1:6] == "ENDBL":
                mode = ""
        elif mode == "pos":
            cols = line.split()
            if cols[0] == "H":
                Hs.append([float(x) for x in cols[1:4]])
            elif cols[0] == "O":
                Os.append([float(x) for x in cols[1:4]])
        elif mode == "lat":
            cols = line.split()
            cell.append([float(x) for x in cols])
    Hs = np.array(Hs)
    Os = np.array(Os)
    cell = np.array(cell)
    return cell, Os, Hs


engel2018 = "Engel, E.A., Anelli, A., Ceriotti, M. et al. Mapping uncharted territory in ice from zeolite networks to ice structures. Nat Commun 9, 2173 (2018). https://doi.org/10.1038/s41467-018-04618-6"


def make_plugin(cell, Os, number, name):
    desc = {"ref": {name: engel2018,
                    "engel{0:02d}".format(number): engel2018},
            "usage": "No options available.",
            "brief": "Hypothetical zeolitic ice"}

    lines = []
    lines.append("desc=" + json.dumps(desc, indent=4))
    lines.append("import numpy as np")
    lines.append("import genice2.lattices")
    lines.append("from genice2.cell import cellvectors")
    lines.append("")
    lines.append("class Lattice(genice2.lattices.Lattice):")
    lines.append("    def __init__(self):")
    lines.append("        self.cell = np.array([")
    for v in cell:
        lines.append("            [{0}, {1}, {2}],".format(*v))
    lines.append("        ])")
    lines.append("        self.waters = np.array([")
    for v in Os:
        lines.append("            [{0}, {1}, {2}],".format(*v))
    lines.append("        ])")
    lines.append("        self.coord = 'relative'")
    lines.append("")
    return "\n".join(lines)


ices = {
    1: "207_1_4435",
    2: "12_2_29187",
    3: "ACO",
    4: "LTA",
    5: "BSV",
    6: "169_2_7915",
    7: "53_3_726600",
    8: "20_2_26425",
    9: "12_2_32449",
    10: "84_2_1419",
    11: "61_2_8842",
    12: "169_2_10608",
    13: "PCOD8047078",
    14: "67_2_1563",  # missing in SI
    15: "PCOD8172143",
    16: "152_2_118474",
    17: "DDR",
    18: "11_2_15848",
    19: "91_2_8335121",
    20: "PCOD8301974",  # missing in SI
    21: "PCOD8045578",
    22: "58_2_511",
    23: "151_2_4949650",  # missing in SI
    24: "PCOD8007225",  # missing in SI
    25: "2_2_342692",
    26: "PCOD8321499",  # missing in SI
    27: "PCOD8047931",
    28: "15_2_201714",
    29: "MAR",
    30: "PCOD8324623",  # missing in SI
    31: "SGT",
    32: "20_2_28176",
    33: "14_2_48453",
    34: "NON",
}


assert False, "Do not run it. Some files are modified by hand."
for num, ice in ices.items():
    try:
        with open(ice + ".castep") as file:
            cell, Os, Hs = load_castep(file)
            outname = "../engel{0:02d}.py".format(num)
            if not os.path.exists(outname):
                print(outname)
                with open(outname, "w") as outfile:
                    outfile.write(make_plugin(cell, Os, num, ice))
    except BaseException:
        if len(ice) == 3:
            try:
                os.symlink(ice + ".py", "../engel{0:02d}.py".format(num))
            except BaseException:
                pass
        pass
