import numpy as np
import genice3.molecule

desc = {"usage": "No options available.", "brief": "An all-atom methane model."}


class Molecule(genice3.molecule.Molecule):
    """
    メタン分子を定義するクラス。
    """

    def __init__(self):
        sites_ = np.array(
            [
                [0.0, 0.0, 0.0],
                [-1.0, -1.0, -1.0],
                [-1.0, +1.0, +1.0],
                [+1.0, -1.0, +1.0],
                [+1.0, +1.0, -1.0],
            ]
        )  # CHHHH
        CH = 0.109  # nm
        sites = sites_ * CH / (3.0**0.5)

        labels = ["C", "H", "H", "H", "H"]
        name = "CH4"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)

