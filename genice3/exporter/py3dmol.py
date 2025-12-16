import py3Dmol
from genice3.plugin import Exporter
import sys
from io import TextIOWrapper
from genice3.genice import GenIce3


def show(genice: GenIce3, **kwargs):

    gro = Exporter("gromacs").dumps(genice, **kwargs)
    view = py3Dmol.view()
    view.addModel(gro, "gro")
    view.setStyle(
        {
            "stick": {},
        }
    )
    view.addUnitCell()
    view.zoomTo()
    return view


def dump(genice: GenIce3, file: TextIOWrapper = sys.stdout, **options):
    raise NotImplementedError("dump is not implemented")
