# coding: utf-8
"""
Show depolarization process of GenIce algorithm in Yaplot format.
defined in https://github.com/vitroid/Yaplot
"""

def hook4(lattice):
    lattice.logger.info("Hook4: Output depolarization process in Yaplot format.")
    print(lattice.yapresult, end="")
    lattice.logger.info("Hook4: end.")
    

hooks = {4:hook4}
