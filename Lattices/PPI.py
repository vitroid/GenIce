import numpy as np

density = 1.7125727202     #default density


bondlen = 3.10      #bond threshold for too much bonds
bondlen = 2.80      #bond threshold for insufficient bonds
celltype = "rect"
cell = """
7.1565 7.1565 7.1565
"""




#Methane A structure; H.E.Maynard-Casely et al, JCP 133, 064504 (2010); doi:10.1063/1.3455889
waters = [[0.186,0.186,0.186],  #C3
          [0.445,0.318,0.839],  #C5
          [0.938,0.072,0.420],  #C8
          [0.839,0.445,0.318],  #C5 [3]
          [0.420,0.938,0.072],  #C8
          [0.318,0.839,0.445],  #C5
          [0.072,0.420,0.938],  #C8
          [0.519,0.189,0.435],  #C7
          [0.204,0.085,0.703],  #C9
          [0.435,0.519,0.189],  #C7
          [0.703,0.204,0.085],  #C9 [10]
          [0.189,0.435,0.519],  #C7
          [0.085,0.703,0.204],  #C9
          [0.575,0.575,0.575],  #C1 [13]
          [0.833,0.291,0.701],  #C4
          [0.701,0.833,0.291],  #C4
          [0.291,0.701,0.833],  #C4
          [0.720,0.626,0.953],  #C6
          [0.953,0.720,0.626],  #C6
          [0.626,0.953,0.720],  #C6
          [0.953,0.954,0.953]]  #C2

waters = np.array(waters)
origin = waters[13].copy()
for water in waters:
    water -= origin
    water -= np.floor(water)


coord = "relative"
