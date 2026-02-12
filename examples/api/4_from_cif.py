from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter
from logging import basicConfig, INFO

# corresponding command:
# genice3 'CIF[file=cif/MEP.cif, osite=T]' --exporter gromacs
basicConfig(level=INFO)
genice = GenIce3()
genice.unitcell = UnitCell("CIF", file="cif/MEP.cif", osite="T")
Exporter("gromacs").dump(genice)
