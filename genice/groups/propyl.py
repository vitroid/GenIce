from genice import group_helper as helper


# This is just an experimental implementation.
def arrange_atoms(cpos, root, cell, molname):
    return helper.Alkyl(cpos, root, cell, molname, ["Ma", ["Mb", "Mc"]])
