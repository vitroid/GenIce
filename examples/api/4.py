from genice3.genice import GenIce3
from genice3.plugin import UnitCell, Exporter

# corresponding command:
# genice3 'CIF[file=cif/MEP.cif, osite=T]' --exporter gromacs
genice = GenIce3()
genice.unitcell = UnitCell("CIF", file="cif/MEP.cif", osite="T")
Exporter("gromacs").dump(genice)
