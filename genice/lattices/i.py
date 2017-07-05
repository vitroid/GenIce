# coding: utf-8
"""
Hypothetical ice "i".

Data source: Fennell, C. J. & Gezelter, J. D. Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice I h. J. Chem. Theory Comput. 1, 662-667 (2005).
"""

density = 0.92     #default density

bondlen = 0.4      #bond threshold	 
celltype = "rect"

cell = """
1.0 1.0 0.94
"""

waters = """
0.16666 0.16666 0.0
0.16666 0.16666 0.47
0.16666 0.83333 0.0
0.16666 0.83333 0.47
0.83333 0.16666 0.0
0.83333 0.16666 0.47
0.83333 0.83333 0.0
0.83333 0.83333 0.47
0.33333 0.33333 0.235
0.33333 0.33333 0.705
0.33333 0.66666 0.235
0.33333 0.66666 0.705
0.66666 0.33333 0.235
0.66666 0.33333 0.705
0.66666 0.66666 0.235
0.66666 0.66666 0.705
"""

coord = "absolute"
