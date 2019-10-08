#!/usr/bin/python
"""

Computational Materials Science 49 (2010) S170–S175
Predicting the hydrogen bond ordered structures of ice Ih, II, III, VI and ice VII: DFT methods with localized based set
Xiaofeng Fan a, Dan Bing a, Jingyun Zhang a, Zexiang Shen a, Jer-Lai Kuo a,b,*

Figure 1(d) #19
"""
import numpy as np

desc={"ref": {"[11#19]": "Xiaofeng Fan, Dan Bing, Jingyun Zhang, Zexiang Shen, Jer-Lai Kuo, Computational Materials Science 49 (2010) S170–S175"},
      "usage": "No options available.",
      "brief": "A candidate for an antiferroelectric Ice XI #19."
      }



waters="""
0.25 0     0.125
0.25 0.333 0
0.25 0     0.5
0.25 0.333 0.625
0.75 0     0.125
0.75 0.333 0
0.75 0     0.5
0.75 0.333 0.625
0    0.5   0.125
0    0.833 0
0    0.5   0.5
0    0.833 0.625
0.5  0.5   0.125
0.5  0.833 0
0.5  0.5   0.5
0.5  0.833 0.625
"""

waters = np.fromstring(waters, dtype=float, sep=' ').reshape(16,3)

fixed="""
0 1
0 13
1 3
1 12
2 11
2 0
3 2
3 10
4 5
4 9
5 8
5 7
6 15
6 4
7 6
7 14
8 1
8 10
9 8
9 0
10 11
10 7
11 6
11 9
12 5
12 14
13 12
13 4
14 15
14 3
15 2
15 13
"""

pairs = fixed

a = 4.4923/10*2
b = 7.7808/10
c = 7.3358/10

from genice.cell import cellvectors
cell  = cellvectors(a,b,c)

density = 18*16/6.022e23 / (np.linalg.det(cell)*1e-21)

coord = "relative"
