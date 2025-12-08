"""
[1] J. L. F. Abascal, E. Sanz, R. G. Fernández, and C. Vega, A potential model for the study of ices and amorphous water: TIP4P/Ice, J. Chem. Phys. 122 (2005) 234511.
"""

from math import pi, sin, cos
import genice3.molecule
import numpy as np


class Molecule(genice3.molecule.Molecule):
    """
    TIP4P/Iceモデルの水分子を定義するクラス。
    """

    def __init__(self):
        L1 = 0.9572 / 10
        L2 = 0.1577 / 10
        theta = 104.52 * pi / 180

        hy = L1 * sin(theta / 2)
        hz = L1 * cos(theta / 2)
        mz = L2
        sites = np.array(
            [[0.0, 0.0, 0.0], [0.0, hy, hz], [0.0, -hy, hz], [0.0, 0.0, mz]]
        )
        sites -= (sites[1] + sites[2] + sites[0] * 0) / 18

        labels = ["OW", "HW1", "HW2", "MW"]
        name = "ICE"
        is_water = True
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)


if __name__ == "__main__":
    water = Molecule()
    print(np.array([16, 1, 1, 0]) @ water.sites)

