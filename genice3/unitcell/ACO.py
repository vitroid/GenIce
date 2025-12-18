"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[ACO] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""

desc = {'ref': {'engel03': 'Engel 2018', 'ACO': 'IZA Database'}, 'usage': 'No options available.', 'brief': 'Hypothetical zeolitic ice'}

import genice3.unitcell
import numpy as np
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    ACO単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.3436    0.3436    0.6564
        0.8436    0.8436    0.1564
        0.3436    0.3436    0.3436
        0.8436    0.8436    0.8436
        0.6564    0.3436    0.3436
        0.1564    0.8436    0.8436
        0.6564    0.3436    0.6564
        0.1564    0.8436    0.1564
        0.3436    0.6564    0.3436
        0.8436    0.1564    0.8436
        0.3436    0.6564    0.6564
        0.8436    0.1564    0.1564
        0.6564    0.6564    0.3436
        0.1564    0.1564    0.8436
        0.6564    0.6564    0.6564
        0.1564    0.1564    0.1564
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 0.30360000000000015

        # density = 0.6961850990811413

        cell = cellvectors(a=0.88235294, b=0.88235294, c=0.88235294)

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )