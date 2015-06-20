import numpy as np
density = 0.92     #default density


bondlen = 1.9      #bond threshold	 
celltype = "rect"
cell = np.array([4.,4.,4.])

positions = np.fromstring("""
0 0 0
2 2 0
2 0 2
0 2 2
1 1 1
3 3 1
3 1 3
1 3 3
""", sep=" ")
positions = positions.reshape((positions.size//3,3))
