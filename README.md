# GenIce
Generate proton-disordered ice structures for GROMACS.

##Requirements
* Python 3
* NetworkX
* numpy

Note: WinPython includes all of these requirements.
##Usage
    usage: genice [-h] [--rep REP REP REP] [--dens DENS] [--seed SEED]
                  [--format gmeqdX] [--water model] [--g12 model] [--g14 model]
                  [--g15 model] [--g16 model] [--cages]
                  Type
    
    positional arguments:
      Type                  Crystal type (1c,1h,etc.)
    
    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [2,2,2]
      --dens DENS, -d DENS  Specify the ice density in g/cm3
      --seed SEED, -s SEED  Random seed [1000]
      --format gmeqdX, -f gmeqdX
                            Specify file format
                            [g(romacs)|m(dview)|e(uler)|q(uaternion)|d(igraph)]
      --water model, -w model
                            Specify water model. (tip3p, tip4p, etc.)
      --g12 model, -D model
                            Specify guest in the 12-hedral cage. (empty, co2,
                            uathf, etc.)
      --g14 model, -T model
                            Specify guest in the 14-hedral cage. (empty, co2,
                            uathf, etc.)
      --g15 model, -P model
                            Specify guest in the 15-hedral cage. (empty, co2,
                            uathf, etc.)
      --g16 model, -H model
                            Specify guest in the 16-hedral cage. (empty, co2,
                            uathf, etc.)
      --cages, -c           Also output the cage positions. (g or m format only)


##Example

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

    ./genice --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format:

    ./genice --cages --g12 co2 --g14 co2 --water tip4p CS1 > cs1.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing
THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

    ./genice --cages --g12 empty --g16 uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro

##Structure generation
The program generates various ice lattice with proton disorder and
without defect.  Total dipole moment is always set to zero.  The
minimal structure (with --rep 1 1 1 option) is not always the unit
cell of the lattice because it is difficult to deal with the hydrogen
bond network topology of tiny lattice under periodic boundary
condition.  Note that the generated structure is not optimal according
to the potential energy.

##Ice structures
* Most popular Ices 1h and Ic.
* Hypothetical Proton-disordered Ice 2.
* Conventional high-pressure ices 3, 4,  6, 7, and 12.
* Monoclinic ice 5 (testing).
* Negative-pressure ice 16.  [Falenty, A., Hansen, T. C. & Kuhs, W. F. Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature 516, 231-233 (2014).]
* Hypothetical ice 0.  [Russo, J., Romano, F. & Tanaka, H. New metastable form of ice and its role in the homogeneous crystallization of water. Nat Mater 13, 733-739 (2014).]
* Hypothetical ice i.  [Fennell, C. J. & Gezelter, J. D. Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice I h. J. Chem. Theory Comput. 1, 662-667 (2005).]
* Hypothetical ices sT' and C0-II.  [Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560-3564 (2013).]
* Clathrate hydrates CS1, CS2, HS1, and TS1.  [Matsumoto, M. & Tanaka, H. On the structure selectivity of clathrate hydrates. J. Phys. Chem. B 115, 8257-8265 (2011).]

Please ask vitroid@gmail.com to add new ice structures.
##Water models
* 3-site: TIP3P (default)
* 4-site: TIP4P
* 5-site: TIP5P

##追記(In preparation)
一部の単位胞(Lattice/1h_unit.pyなど)は、その大きさが小さすぎて、グラフを定義できないため、単位胞として2x1x1倍格子(1h.py)を収録している。しかし、2x1x1単位胞をrepeatすると、x軸方向が単位格子の奇数倍の格子を作れない。

そのような場合のために、単位胞の座標だけを定数倍した、新しい単位胞をpython moduleの形で作る機能を追加した。例えば、
    genice --format X --density 0.92 -r 7 1 1 1h_unit > 1hx711
で、単位胞の7x1x1倍の構造をpython module形式で生成できる。これを使って、グラフを含む7x4x5倍格子を作りたい場合は、通常通り、
    genice --format g -r 1 4 5 1hx711 > 1hx745.gro
などとすれば良い。
