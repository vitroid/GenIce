"""

Command line: /Users/matto/anaconda3/envs/genice/bin/genice2 5 -r 1 1 3 -f reshape[1,0,0,0,1,0,1,0,1]
Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[1 0 1]
"""

from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {},
        "usage": "No options available.",
        "brief": "Ice V with orthogonal unit cell. (testing)"
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 3.00000000000215
        self.coord = 'relative'
        self.density = 1.23983176668
        self.waters = """
            0.7375    0.1404    0.1618
            0.4292    0.1404    0.1715
            0.9292    0.8596    0.1715
            0.2375    0.8596    0.1618
            0.5708    0.6404    0.3285
            0.5958    0.6404    0.0049
            0.0958    0.3596    0.0049
            0.0708    0.3596    0.3285
            0.5259    0.8475    0.2492
            0.6408    0.8475    0.0841
            0.1408    0.1525    0.0841
            0.0259    0.1525    0.2492
            0.6925    0.3475    0.0826
            0.4741    0.3475    0.2508
            0.9741    0.6525    0.2508
            0.1925    0.6525    0.0826
            0.7448    0.4435    0.2181
            0.4219    0.4435    0.1152
            0.9219    0.5565    0.1152
            0.2448    0.5565    0.2181
            0.9114    0.9435    0.0515
            0.2552    0.9435    0.2819
            0.7552    0.0565    0.2819
            0.4114    0.0565    0.0515
            0.5833    0.6847    0.1667
            0.0833    0.3153    0.1667
            0.7500    0.1847    0.0000
            0.2500    0.8153    0.0000
            0.4042    0.1404    0.4951
            0.0958    0.1404    0.5049
            0.5958    0.8596    0.5049
            0.9042    0.8596    0.4951
            0.2375    0.6404    0.6618
            0.2625    0.6404    0.3382
            0.7625    0.3596    0.3382
            0.7375    0.3596    0.6618
            0.1925    0.8475    0.5826
            0.3075    0.8475    0.4174
            0.8075    0.1525    0.4174
            0.6925    0.1525    0.5826
            0.3592    0.3475    0.4159
            0.1408    0.3475    0.5841
            0.6408    0.6525    0.5841
            0.8592    0.6525    0.4159
            0.4114    0.4435    0.5515
            0.0886    0.4435    0.4485
            0.5886    0.5565    0.4485
            0.9114    0.5565    0.5515
            0.5781    0.9435    0.3848
            0.9219    0.9435    0.6152
            0.4219    0.0565    0.6152
            0.0781    0.0565    0.3848
            0.2500    0.6847    0.5000
            0.7500    0.3153    0.5000
            0.4167    0.1847    0.3333
            0.9167    0.8153    0.3333
            0.0708    0.1404    0.8285
            0.7625    0.1404    0.8382
            0.2625    0.8596    0.8382
            0.5708    0.8596    0.8285
            0.9042    0.6404    0.9951
            0.9292    0.6404    0.6715
            0.4292    0.3596    0.6715
            0.4042    0.3596    0.9951
            0.8592    0.8475    0.9159
            0.9741    0.8475    0.7508
            0.4741    0.1525    0.7508
            0.3592    0.1525    0.9159
            0.0259    0.3475    0.7492
            0.8075    0.3475    0.9174
            0.3075    0.6525    0.9174
            0.5259    0.6525    0.7492
            0.0781    0.4435    0.8848
            0.7552    0.4435    0.7819
            0.2552    0.5565    0.7819
            0.5781    0.5565    0.8848
            0.2448    0.9435    0.7181
            0.5886    0.9435    0.9485
            0.0886    0.0565    0.9485
            0.7448    0.0565    0.7181
            0.9167    0.6847    0.8333
            0.4167    0.3153    0.8333
            0.0833    0.1847    0.6667
            0.5833    0.8153    0.6667
        """

        self.cell = cellvectors(a=9.19977813,
                                b=7.52346281,
                                c=29.25857353)
