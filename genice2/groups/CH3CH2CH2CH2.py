import genice2.groups
from genice2.alkyl import Alkyl


class Group(genice2.groups.Group):
    def arrange_atoms(self, cage_center, root_position, cell, molname, origin_atom):
        """
        put a butyl group rooted at root_position toward cage_center.
        """
        bondlen={("C", "C"): 0.153,
                ("C", "N"): 0.1471,
                ("N", "C"): 0.1471,
                ("H", "C"): 0.1090,
                ("C", "H"): 0.1090,
                }
        return Alkyl(cage_center, root_position, cell, molname, ["C1", ["C2", ["C3", ["C4", "H", "H", "H"], "H", "H"], "H", "H"], "H", "H"], bondlen, origin_atom=origin_atom)
