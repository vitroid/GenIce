# coding: utf-8
"""
Note: Due to the technical limitation in the GenIce algorithm, the minimum lattice size is larger than the crystallographic unit cell size.
"""

density = 1.6     #default density

celltype = "rect"
cell = """
4 4 4
"""

waters = """
0.00 0.00 0.00
0.50 0.50 0.00
0.50 0.00 0.50
0.00 0.50 0.50
0.25 0.25 0.25
0.75 0.75 0.25
0.75 0.25 0.75
0.25 0.75 0.75
0.50 0.00 0.00
0.00 0.50 0.00
0.00 0.00 0.50
0.50 0.50 0.50
0.75 0.25 0.25
0.25 0.75 0.25
0.25 0.25 0.75
0.75 0.75 0.75
"""
coord = "relative"
double_network = True  #It is necessary only for ices 6 and 7

pairs = """
0 4
0 6
1 4
1 5
1 6
2 6
3 4
3 5
3 6
3 7
4 2
5 0
5 2
7 0
7 1
7 2
8 12
8 14
9 12
9 13
9 14
10 14
11 12
11 13
11 14
11 15
12 10
13 8
13 10
15 8
15 9
15 10
"""

