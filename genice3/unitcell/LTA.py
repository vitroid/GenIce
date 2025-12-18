"""

Command line: /Volumes/workarea/venvs/genice2/bin/genice zeolite[LTA] -f python
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""

desc = {'ref': {'engel04': 'Engel 2018', 'LTA': 'IZA Database'}, 'usage': 'No options available.', 'brief': 'Hypothetical zeolitic ice'}

import genice3.unitcell
import numpy as np
from cif2ice import cellvectors


class UnitCell(genice3.unitcell.UnitCell):
    """
    LTA単位胞を定義するクラス。
    """

    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成

        waters = np.fromstring(
            """
        0.6316    0.1823    0.0000
        0.8177    0.0000    0.3684
        0.6316    0.8177    0.0000
        0.8177    0.0000    0.6316
        0.1823    0.6316    0.0000
        0.0000    0.8177    0.3684
        0.8177    0.6316    0.0000
        0.0000    0.8177    0.6316
        0.0000    0.6316    0.8177
        0.3684    0.8177    0.0000
        0.0000    0.6316    0.1823
        0.0000    0.3684    0.8177
        0.1823    0.0000    0.6316
        0.3684    0.1823    0.0000
        0.0000    0.3684    0.1823
        0.1823    0.0000    0.3684
        0.6316    0.0000    0.8177
        0.8177    0.3684    0.0000
        0.6316    0.0000    0.1823
        0.3684    0.0000    0.8177
        0.0000    0.1823    0.6316
        0.1823    0.3684    0.0000
        0.3684    0.0000    0.1823
        0.0000    0.1823    0.3684
        """,
            sep=" ",
        ).reshape(-1, 3)

        coord = "relative"

        bondlen = 0.30360000000000004

        # density = 0.5846833799862865

        cell = cellvectors(a=1.07055113, b=1.07055113, c=1.07055113)

        super().__init__(
            cell=cell,
            lattice_sites=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )