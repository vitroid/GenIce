# coding: utf-8
"""
Data source:

Petrenko and Whitworth, Physics of Ice, Table 11.2
"""
density = 1.373     #default density


bondlen = 3      #bond threshold	 
celltype = "rect"
cell = """
6.181 6.181 5.698
"""

double_network = True  #It is necessary only for ices 6 and 7

coord = "relative"

waters = """
    0.7504    0.2502    0.7510
    0.4707    0.2502    0.3666
    0.9710    0.7509    0.6348
    0.5298    0.7509    0.6348
    0.2501    0.4710    0.8674
    0.0295    0.2502    0.3666
    0.7504    0.5301    0.1341
    0.2501    0.7509    0.2503
    0.7504    0.9716    0.1341
    0.2501    0.0295    0.8674
"""

pairs = """
1 0
3 2
4 2
4 3
5 0
5 1
6 0
6 1
6 5
7 2
7 3
7 4
8 0
8 1
8 5
8 6
9 2
9 3
9 4
9 7
"""
