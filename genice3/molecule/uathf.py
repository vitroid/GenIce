import numpy as np

# United-atom THF model

import genice3.molecule


desc = {
    "usage": "No options available.",
    "brief": "A united-atom five-site tetrahydrofuran (THF) model.",
}


class Molecule(genice3.molecule.Molecule):
    """
    統合原子5サイトテトラヒドロフラン（THF）モデルを定義するクラス。
    """

    def __init__(self):
        sites = np.array(
            [
                [0.0, -0.119625, 0.0],
                [0.116284, -0.039705, 0.0],
                [0.076453, 0.107915, 0.0],
                [-0.076453, 0.107915, 0.0],
                [-0.116284, -0.039705, 0.0],
            ]
        )

        labels = ["O", "CA", "CB", "CB", "CA"]
        name = "THF"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)

