from genice import group_helper as helper


# This is just an experimental implementation.
def arrange_atoms(cpos, root, cell, molname):
    """
    2,2-dimethylpropyl group rooted at root toward cpos.
    It interferes with host water molecules.
    """
    return helper.Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", "Mc", "Md", "Me"]])

