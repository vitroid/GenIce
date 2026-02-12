# 分極した氷の作り方 (1) コンストラクタで target_pol を指定
# corresponding command: 7_polarized_1.sh または 7_polarized_1_flat.sh

from logging import basicConfig, INFO
import numpy as np
from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter

basicConfig(level=INFO)

genice = GenIce3(
    seed=114,
    depol_loop=1000,
    replication_matrix=np.diag([2, 2, 2]),
    target_pol=np.array([4.0, 0.0, 0.0]),
)
genice.unitcell = UnitCell("1h")

Exporter("_pol").dump(genice)
