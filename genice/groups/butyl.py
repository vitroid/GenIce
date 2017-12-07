"""
put a butyl group rooted at root toward cpos.
"""
from genice import group_helper as helper


# This is just an experimental implementation.
def arrange_atoms(cpos, root, cell, molname):
    tree = ["Ma", ["Mb", ["Mc", "Md"]]]
    return helper.Alkyl(cpos, root, cell, molname, tree)