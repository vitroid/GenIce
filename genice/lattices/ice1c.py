# coding: utf-8
"""
Data sources
[C2] W. L. Vos, L. W. Finger, R. J. Hemley, H.-K. Mao, Phys. Rev. Lett., 1993, DOI:10.1103/PhysRevLett.71.3150.
"""
density = 0.92     #default density


bondlen = 1.9      #bond threshold	 
celltype = "rect"
cell = """
4 4 4
"""

waters = """
0 0 0
0.5 0.5 0
0.5 0 0.5
0 0.5 0.5
0.25 0.25 0.25
0.75 0.75 0.25
0.75 0.25 0.75
0.25 0.75 0.75
"""
coord = "relative"

pairs="""
0 4
0 5
0 6
0 7
1 4
1 5
1 6
1 7
2 4
2 5
2 6
2 7
3 4
3 5
3 6
3 7
"""
