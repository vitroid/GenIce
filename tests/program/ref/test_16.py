"""

Data source:
[ice 16] Falenty, A., Hansen, T. C. & Kuhs, W. F. Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature 516, 231-233 (2014).
[C15] Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.
[sII] Jeffrey, G A. “Hydrate Inclusion Compounds.” Inclusion Compounds 1 (1984): 135–190.
[CS2] Kosyakov, Viktor I, and T M Polyanskaya. “Using Structural Data for Estimating the Stability of Water Networks in Clathrate and Semiclathrate Hydrates.” Journal of Structural Chemistry 40.2 (1999): 239–245.
[MTN] http://www.iza-structure.org/databases/

Reshaping the unit cell.
  i:[1 0 0]
  j:[0 1 0]
  k:[0 0 1]
"""
bondlen=3.034768428692841
coord='relative'
celltype='rect'
cell='17.12097716 17.12097716 34.24195431'
density=0.81
waters="""
    0.5680    0.8708    0.2840
    0.6820    0.3708    0.0910
    0.3180    0.9320    0.0604
    0.9320    0.6292    0.2160
    0.9320    0.9320    0.0646
    0.3708    0.6820    0.0910
    0.3180    0.6208    0.2160
    0.8180    0.4320    0.0604
    0.8750    0.3750    0.4375
    0.3180    0.6292    0.4090
    0.4320    0.1208    0.4090
    0.2180    0.2180    0.3590
    0.6292    0.4320    0.4660
    0.6208    0.4320    0.1590
    0.1292    0.8180    0.4090
    0.2820    0.7820    0.3910
    0.4320    0.4320    0.0646
    0.2180    0.0320    0.2660
    0.6208    0.8180    0.4660
    0.1820    0.8792    0.2840
    0.1820    0.0680    0.1896
    0.4680    0.4680    0.1410
    0.9320    0.1208    0.1590
    0.8750    0.8750    0.1875
    0.1292    0.4320    0.2160
    0.1820    0.8708    0.0910
    0.4320    0.9320    0.3146
    0.1208    0.3180    0.4660
    0.6820    0.0680    0.4396
    0.8180    0.8180    0.0646
    0.8792    0.0680    0.3410
    0.9680    0.9680    0.1410
    0.0680    0.6820    0.4396
    0.7180    0.7180    0.3590
    0.8708    0.1820    0.0910
    0.5320    0.5320    0.3590
    0.6820    0.8708    0.3410
    0.0680    0.8792    0.3410
    0.5680    0.3792    0.3410
    0.3180    0.3180    0.0646
    0.3180    0.1292    0.1590
    0.1208    0.4320    0.4090
    0.8708    0.0680    0.0340
    0.0680    0.0680    0.4354
    0.8792    0.1820    0.2840
    0.2820    0.2820    0.1410
    0.3708    0.5680    0.0340
    0.1250    0.6250    0.0625
    0.8708    0.6820    0.3410
    0.9680    0.4680    0.3910
    0.6820    0.6820    0.4354
    0.6820    0.3792    0.2840
    0.4680    0.7820    0.4840
    0.0680    0.1820    0.1896
    0.3180    0.4320    0.3104
    0.4320    0.6292    0.4660
    0.5320    0.2180    0.0160
    0.6292    0.9320    0.2160
    0.4320    0.8180    0.0604
    0.7180    0.2180    0.1090
    0.5680    0.3708    0.0340
    0.1820    0.6820    0.1854
    0.6820    0.1820    0.1854
    0.0680    0.3792    0.0910
    0.9320    0.3180    0.0604
    0.0680    0.5680    0.1854
    0.9680    0.2820    0.4840
    0.5320    0.0320    0.1090
    0.8180    0.6292    0.1590
    0.3708    0.0680    0.2840
    0.1292    0.3180    0.1590
    0.8708    0.5680    0.2840
    0.9320    0.1292    0.4660
    0.3750    0.3750    0.1875
    0.2820    0.9680    0.4840
    0.1820    0.3708    0.3410
    0.8180    0.9320    0.3104
    0.3750    0.8750    0.4375
    0.8180    0.6208    0.4660
    0.6250    0.1250    0.0625
    0.6250    0.6250    0.3125
    0.1820    0.5680    0.4396
    0.5680    0.6820    0.1896
    0.7820    0.9680    0.2340
    0.1208    0.9320    0.1590
    0.1820    0.1820    0.4354
    0.3792    0.0680    0.0910
    0.8792    0.6820    0.0340
    0.7820    0.4680    0.4840
    0.7180    0.5320    0.2660
    0.9320    0.6208    0.4090
    0.0680    0.3708    0.2840
    0.2180    0.7180    0.1090
    0.6820    0.8792    0.0340
    0.8180    0.1208    0.2160
    0.9680    0.7820    0.2340
    0.6292    0.3180    0.4090
    0.4320    0.3180    0.3104
    0.3180    0.1208    0.4660
    0.4680    0.2820    0.2340
    0.1292    0.9320    0.4660
    0.0320    0.7180    0.0160
    0.0680    0.8708    0.0340
    0.4680    0.9680    0.3910
    0.0320    0.0320    0.3590
    0.9320    0.4320    0.3146
    0.4320    0.1292    0.2160
    0.9320    0.8180    0.3104
    0.2820    0.4680    0.2340
    0.6208    0.9320    0.4090
    0.3708    0.1820    0.3410
    0.5320    0.7180    0.2660
    0.3792    0.5680    0.3410
    0.3792    0.1820    0.0340
    0.8180    0.1292    0.4090
    0.7180    0.0320    0.0160
    0.2180    0.5320    0.0160
    0.3792    0.6820    0.2840
    0.5680    0.1820    0.4396
    0.0320    0.2180    0.2660
    0.7820    0.7820    0.1410
    0.5680    0.0680    0.1854
    0.4320    0.6208    0.1590
    0.5680    0.8792    0.0910
    0.6292    0.8180    0.1590
    0.5680    0.5680    0.4354
    0.8792    0.5680    0.0910
    0.0320    0.5320    0.1090
    0.8180    0.3180    0.3146
    0.6208    0.3180    0.2160
    0.1250    0.1250    0.3125
    0.1820    0.3792    0.0340
    0.3180    0.8180    0.3146
    0.1208    0.8180    0.2160
    0.7820    0.2820    0.3910
    0.6820    0.5680    0.1896
    0.5680    0.8708    0.7840
    0.6820    0.3708    0.5910
    0.3180    0.9320    0.5604
    0.9320    0.6292    0.7160
    0.9320    0.9320    0.5646
    0.3708    0.6820    0.5910
    0.3180    0.6208    0.7160
    0.8180    0.4320    0.5604
    0.8750    0.3750    0.9375
    0.3180    0.6292    0.9090
    0.4320    0.1208    0.9090
    0.2180    0.2180    0.8590
    0.6292    0.4320    0.9660
    0.6208    0.4320    0.6590
    0.1292    0.8180    0.9090
    0.2820    0.7820    0.8910
    0.4320    0.4320    0.5646
    0.2180    0.0320    0.7660
    0.6208    0.8180    0.9660
    0.1820    0.8792    0.7840
    0.1820    0.0680    0.6896
    0.4680    0.4680    0.6410
    0.9320    0.1208    0.6590
    0.8750    0.8750    0.6875
    0.1292    0.4320    0.7160
    0.1820    0.8708    0.5910
    0.4320    0.9320    0.8146
    0.1208    0.3180    0.9660
    0.6820    0.0680    0.9396
    0.8180    0.8180    0.5646
    0.8792    0.0680    0.8410
    0.9680    0.9680    0.6410
    0.0680    0.6820    0.9396
    0.7180    0.7180    0.8590
    0.8708    0.1820    0.5910
    0.5320    0.5320    0.8590
    0.6820    0.8708    0.8410
    0.0680    0.8792    0.8410
    0.5680    0.3792    0.8410
    0.3180    0.3180    0.5646
    0.3180    0.1292    0.6590
    0.1208    0.4320    0.9090
    0.8708    0.0680    0.5340
    0.0680    0.0680    0.9354
    0.8792    0.1820    0.7840
    0.2820    0.2820    0.6410
    0.3708    0.5680    0.5340
    0.1250    0.6250    0.5625
    0.8708    0.6820    0.8410
    0.9680    0.4680    0.8910
    0.6820    0.6820    0.9354
    0.6820    0.3792    0.7840
    0.4680    0.7820    0.9840
    0.0680    0.1820    0.6896
    0.3180    0.4320    0.8104
    0.4320    0.6292    0.9660
    0.5320    0.2180    0.5160
    0.6292    0.9320    0.7160
    0.4320    0.8180    0.5604
    0.7180    0.2180    0.6090
    0.5680    0.3708    0.5340
    0.1820    0.6820    0.6854
    0.6820    0.1820    0.6854
    0.0680    0.3792    0.5910
    0.9320    0.3180    0.5604
    0.0680    0.5680    0.6854
    0.9680    0.2820    0.9840
    0.5320    0.0320    0.6090
    0.8180    0.6292    0.6590
    0.3708    0.0680    0.7840
    0.1292    0.3180    0.6590
    0.8708    0.5680    0.7840
    0.9320    0.1292    0.9660
    0.3750    0.3750    0.6875
    0.2820    0.9680    0.9840
    0.1820    0.3708    0.8410
    0.8180    0.9320    0.8104
    0.3750    0.8750    0.9375
    0.8180    0.6208    0.9660
    0.6250    0.1250    0.5625
    0.6250    0.6250    0.8125
    0.1820    0.5680    0.9396
    0.5680    0.6820    0.6896
    0.7820    0.9680    0.7340
    0.1208    0.9320    0.6590
    0.1820    0.1820    0.9354
    0.3792    0.0680    0.5910
    0.8792    0.6820    0.5340
    0.7820    0.4680    0.9840
    0.7180    0.5320    0.7660
    0.9320    0.6208    0.9090
    0.0680    0.3708    0.7840
    0.2180    0.7180    0.6090
    0.6820    0.8792    0.5340
    0.8180    0.1208    0.7160
    0.9680    0.7820    0.7340
    0.6292    0.3180    0.9090
    0.4320    0.3180    0.8104
    0.3180    0.1208    0.9660
    0.4680    0.2820    0.7340
    0.1292    0.9320    0.9660
    0.0320    0.7180    0.5160
    0.0680    0.8708    0.5340
    0.4680    0.9680    0.8910
    0.0320    0.0320    0.8590
    0.9320    0.4320    0.8146
    0.4320    0.1292    0.7160
    0.9320    0.8180    0.8104
    0.2820    0.4680    0.7340
    0.6208    0.9320    0.9090
    0.3708    0.1820    0.8410
    0.5320    0.7180    0.7660
    0.3792    0.5680    0.8410
    0.3792    0.1820    0.5340
    0.8180    0.1292    0.9090
    0.7180    0.0320    0.5160
    0.2180    0.5320    0.5160
    0.3792    0.6820    0.7840
    0.5680    0.1820    0.9396
    0.0320    0.2180    0.7660
    0.7820    0.7820    0.6410
    0.5680    0.0680    0.6854
    0.4320    0.6208    0.6590
    0.5680    0.8792    0.5910
    0.6292    0.8180    0.6590
    0.5680    0.5680    0.9354
    0.8792    0.5680    0.5910
    0.0320    0.5320    0.6090
    0.8180    0.3180    0.8146
    0.6208    0.3180    0.7160
    0.1250    0.1250    0.8125
    0.1820    0.3792    0.5340
    0.3180    0.8180    0.8146
    0.1208    0.8180    0.7160
    0.7820    0.2820    0.8910
    0.6820    0.5680    0.6896
"""

cages="""
12    0.5000    0.2500    0.1250
12    0.5000    0.5000    0.2500
12    0.2500    0.5000    0.1250
12    0.0000    0.2500    0.3750
12    0.2500    0.0000    0.3750
12    0.0000    0.5000    0.0000
12    0.2500    0.2500    0.2500
12    0.5000    0.0000    0.0000
12    0.7500    0.5000    0.3750
12    0.0000    0.0000    0.2500
12    0.0000    0.7500    0.1250
12    0.2500    0.7500    0.0000
12    0.7500    0.2500    0.0000
12    0.7500    0.0000    0.1250
12    0.5000    0.7500    0.3750
12    0.7500    0.7500    0.2500
16    0.6250    0.1250    0.3125
16    0.1250    0.6250    0.3125
16    0.8750    0.3750    0.1875
16    0.3750    0.8750    0.1875
16    0.3750    0.3750    0.4375
16    0.6250    0.6250    0.0625
16    0.8750    0.8750    0.4375
16    0.1250    0.1250    0.0625
12    0.5000    0.2500    0.6250
12    0.5000    0.5000    0.7500
12    0.2500    0.5000    0.6250
12    0.0000    0.2500    0.8750
12    0.2500    0.0000    0.8750
12    0.0000    0.5000    0.5000
12    0.2500    0.2500    0.7500
12    0.5000    0.0000    0.5000
12    0.7500    0.5000    0.8750
12    0.0000    0.0000    0.7500
12    0.0000    0.7500    0.6250
12    0.2500    0.7500    0.5000
12    0.7500    0.2500    0.5000
12    0.7500    0.0000    0.6250
12    0.5000    0.7500    0.8750
12    0.7500    0.7500    0.7500
16    0.6250    0.1250    0.8125
16    0.1250    0.6250    0.8125
16    0.8750    0.3750    0.6875
16    0.3750    0.8750    0.6875
16    0.3750    0.3750    0.9375
16    0.6250    0.6250    0.5625
16    0.8750    0.8750    0.9375
16    0.1250    0.1250    0.5625
"""

