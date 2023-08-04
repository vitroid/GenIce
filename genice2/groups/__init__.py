# coding: utf-8
class Group():
    """
    Base class of groups
    """
    def arrange_atoms(self, cage_center, root_position, cell, molname, origin_atom):
        """Place a chemical group rooted at an atom, in a specified cage.

        Args:
            cage_center (numpy rank 3 vector): Center of the cage, in the fractional coordinate
            root_position (numpy rank 3 vector): Position of the root atom, in the fractional coordinate
            cell (numpy array 3x3): Cell matrix.
            molname (str): Name of the chemical group, given as the residue name in Gromacs file.
            origin_atom (str): Name of the root atom. It is required to assume the bond length.
        """
        pass
