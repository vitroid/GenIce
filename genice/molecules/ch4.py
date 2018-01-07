# coding: utf-8
import numpy as np

sites = np.array([[0.0, 0.0, 0.0],
                  [-1.,-1.,-1.],
                  [-1.,+1.,+1.],
                  [+1.,-1.,+1.],
                  [+1.,+1.,-1.]]) #CHHHH
CH = 0.109  # nm
sites *= CH / (3.0**0.5)

labels = ["C","H","H","H","H"]
name = "CH4"
