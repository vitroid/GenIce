

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"III": 'Petrenko 1999'},
        "usage": "No options available.",
        "brief": "Ice III."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 3.0
        self.coord = 'relative'
        self.cell = '6.666 6.666 6.936'
        self.density = 1.165
        self.waters = """
            0.3936    0.3939    0.0141
            0.6103    0.1988    0.4695
            0.8912    0.6997    0.7988
            0.1986    0.6109    0.5452
            0.6991    0.8921    0.2159
            0.6074    0.6080    0.5073
            0.1069    0.8949    0.2536
            0.3019    0.1099    0.7232
            0.3907    0.8032    0.9768
            0.8941    0.1070    0.7609
            0.1098    0.3022    0.2915
            0.8024    0.3911    0.0379
        """

        self.cell = cellvectors(a=6.666,
                                b=6.666,
                                c=6.936)
