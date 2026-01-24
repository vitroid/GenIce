from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter
from logging import basicConfig, INFO

# corresponding command:
# genice3 A15 -e cif
basicConfig(level=INFO)
genice = GenIce3()
genice.unitcell = UnitCell("A15")
Exporter("cif").dump(genice)
