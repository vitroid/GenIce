from logging import basicConfig, DEBUG
import numpy as np
from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter, Molecule

# corresponding command:
# genice "A15[shift=(0.1,0.1,0.1), anion.0=Cl, cation.6=Na, density=0.8]" --rep 2 2 2
# --exporter "gromacs[guest.A12=me, guest.A14=et, spot_guest.0=4site, water=4site]" -D
basicConfig(level=DEBUG)  # -D
genice = GenIce3(replication_matrix=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
me = Molecule("me")
et = Molecule("et")
water = Molecule("4site")
genice.unitcell = UnitCell(
    "A15", shift=(0.1, 0.1, 0.1), anions={0: "Cl"}, cations={6: "Na"}, density=0.8
)
Exporter("gromacs").dump(
    genice, guest={"A12": "me", "A14": "et"}, spot_guest={0: "4site"}, water=water
)
