from genice import group_helper as helper


# This is just an experimental implementation.
def arrange_atoms(cpos, root, cell, molname):
    """
    put a group rooted at root toward cpos.
    """
    tree = ["Ma", ["Mb", ["Mc", "Md", "Me"]]]
    return helper.Alkyl(cpos, root, cell, molname, tree)
