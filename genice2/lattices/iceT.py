from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"T": 'Hirata 2017'},
        "usage": "No options available.",
        "brief": "Hypothetical ice T."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = """
        13.3672639407 13.3672639407 7.42849456948
        """
        self.waters = """
        0.142 0.751 0.622
        0.858 0.249 0.622
        0.251 0.858 0.372
        0.749 0.142 0.372
        0.358 0.749 0.122
        0.642 0.251 0.122
        0.249 0.642 0.872
        0.751 0.358 0.872
        0.327 0.173 0.249
        0.673 0.827 0.249
        0.673 0.673 0.999
        0.327 0.327 0.999
        0.173 0.327 0.749
        0.827 0.673 0.749
        0.827 0.827 0.499
        0.173 0.173 0.499
        0.358 0.251 0.622
        0.642 0.749 0.622
        0.751 0.642 0.372
        0.249 0.358 0.372
        0.142 0.249 0.122
        0.858 0.751 0.122
        0.749 0.858 0.872
        0.251 0.142 0.872
        0.468 0.166 0.971
        0.532 0.834 0.971
        0.666 0.532 0.721
        0.334 0.468 0.721
        0.032 0.334 0.471
        0.968 0.666 0.471
        0.834 0.968 0.221
        0.166 0.032 0.221
        0.532 0.166 0.471
        0.468 0.834 0.471
        0.666 0.468 0.221
        0.334 0.532 0.221
        0.968 0.334 0.971
        0.032 0.666 0.971
        0.834 0.032 0.721
        0.166 0.968 0.721
        0.468 0.334 0.276
        0.532 0.666 0.276
        0.834 0.532 0.026
        0.166 0.468 0.026
        0.032 0.166 0.776
        0.968 0.834 0.776
        0.666 0.968 0.526
        0.334 0.032 0.526
        0.532 0.334 0.776
        0.468 0.666 0.776
        0.834 0.468 0.526
        0.166 0.532 0.526
        0.968 0.166 0.276
        0.032 0.834 0.276
        0.666 0.032 0.026
        0.334 0.968 0.026
        0.5 0.0 0.25
        0.5 0.5 0.0
        0.0 0.5 0.75
        0.0 0.0 0.5
        0.173 0.827 0.999
        0.827 0.173 0.999
        0.327 0.827 0.749
        0.673 0.173 0.749
        0.327 0.673 0.499
        0.673 0.327 0.499
        0.173 0.673 0.249
        0.827 0.327 0.249
        0.5 0.0 0.75
        0.5 0.5 0.5
        0.0 0.5 0.25
        0.0 0.0 0.0
        """
        self.coord = "relative"
        self.bondlen = 2.83
        self.density = 1.62135603834

        self.cell = cellvectors(a=13.3672639407,
                                b=13.3672639407,
                                c=7.42849456948)
