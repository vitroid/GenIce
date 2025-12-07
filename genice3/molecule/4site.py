from math import cos, radians, sin
import numpy as np
import genice3.molecule

desc = {
    "ref": {"TIP4P(a)": "Jorgensen 1983", "TIP4P(b)": "Jorgensen 1985"},
    "usage": "No options available.",
    "brief": "A typical 4-site model.",
}


class Molecule(genice3.molecule.Molecule):
    """
    4サイトモデルの水分子を定義するクラス。
    """

    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.15 / 10
        theta = radians(104.52)

        hy = L1 * sin(theta / 2)
        hz = L1 * cos(theta / 2)
        mz = L2
        sites = np.array(
            [[0.0, 0.0, 0.0], [0.0, hy, hz], [0.0, -hy, hz], [0.0, 0.0, mz]]
        )
        sites -= (sites[1] + sites[2] + sites[3] * 0) / 18
        labels = ["OW", "HW1", "HW2", "MW"]
        name = "ICE"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)
