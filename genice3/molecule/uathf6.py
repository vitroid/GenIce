import numpy as np

# United-atom THF model with a dummy site

import genice3.molecule


class Molecule(genice3.molecule.Molecule):
    """
    ダミーサイトを持つ統合原子THFモデルを定義するクラス。
    """

    def __init__(self):
        sites = np.array(
            [
                [0.0, -0.119625, 0.0],
                [0.116284, -0.039705, 0.0],
                [0.076453, 0.107915, 0.0],
                [-0.076453, 0.107915, 0.0],
                [-0.116284, -0.039705, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )

        # self.atoms_ = ["O", "C", "C", "C", "C", "."]
        labels = ["O", "CA", "CB", "CB", "CA", "CM"]
        name = "THF"
        is_water = False
        super().__init__(sites=sites, labels=labels, name=name, is_water=is_water)

