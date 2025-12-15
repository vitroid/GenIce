import sys

# confirmed:  "8", "iceA", "iceB", "9", "11"
for hydrogen_ordered_ice in ["2"] + sys.argv[1:]:

    from genice2.genice import GenIce
    from genice2.plugin import Lattice, Format, Molecule

    lattice = Lattice(hydrogen_ordered_ice)
    formatter = Format("raw", stage=(1, 2, 3, 4))
    water = Molecule("spce")
    genice2 = GenIce(lattice).generate_ice(formatter, water=water)
    print(genice2["digraph"])

    from genice3.plugin import UnitCell
    from genice3.genice import GenIce3

    genice3 = GenIce3()
    genice3.unitcell = UnitCell(hydrogen_ordered_ice)
    print(genice3.digraph)

    import networkx as nx

    print(f"{hydrogen_ordered_ice}")
    print(nx.is_isomorphic(genice2["graph"], genice3.graph))
    print(nx.is_isomorphic(genice2["digraph"], genice3.digraph))
