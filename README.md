# GenIce
Generate proton-disordered ice structures for GROMACS.

##Requirements

* Python 3
* NetworkX
* numpy

Note: WinPython includes all of these requirements.
##Installation
###Using steuptools (Unix)

1. First, install the [setuptools](https://setuptools.readthedocs.io/en/latest/) by some means.  Skip this step if you already have `pip3` command.
2. Download the source code from [this page](https://github.com/vitroid/GenIce).
3. Run the following command.

    % cd GenIce
    % ./setup.py install

4. Or, just run genice.x command there.

    % cd GenIce
    % ./genice.x 7 > 7.gro 

##Uninstallation

    % pip3 uninstall GenIce
    
##Usage
    usage: genice.x [-h] [--rep REP REP REP] [--dens DENS] [--seed SEED]
                    [--format gmeqdXoc] [--water model] [--guest D=empty]
                    [--debug] [--quiet]
                    Type
    
    positional arguments:
      Type                  Crystal type (1c,1h,etc.)
    
    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [2,2,2]
      --dens DENS, -d DENS  Specify the ice density in g/cm3
      --seed SEED, -s SEED  Random seed [1000]
      --format gmeqdXoc, -f gmeqdXoc
                            Specify file format [g(romacs)|m(dview)|e(uler)|q(uate
                            rnion)|d(igraph)|o(penScad)|c(entersofmass)]
      --water model, -w model
                            Specify water model. (tip3p, tip4p, etc.)
      --guest D=empty, -g D=empty
                            Specify guest in the cage. (D=empty, T=co2, etc.)
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.


##Example

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

    ./genice --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format:

    ./genice -g 12=co2 -g 14=co2 --water tip4p CS1 > cs1.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing
THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

    ./genice -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro

##Structure generation
The program generates various ice lattice with proton disorder and
without defect.  Total dipole moment is always set to zero.  The
minimal structure (with --rep 1 1 1 option) is not always the unit
cell of the lattice because it is difficult to deal with the hydrogen
bond network topology of tiny lattice under periodic boundary
condition.  Note that the generated structure is not optimal according
to the potential energy.

##Ice structures

Symbol | Description| Reference
-------|------------|----------
1h, 1c | Most popular Ice I (hexagonal or cubic)
2d     | Hypothetical Proton-disordered Ice II. |Nakamura, Tatsuya et al. “Thermodynamic Stability of Ice II and Its Hydrogen-Disordered Counterpart: Role of Zero-Point Energy.” The Journal of Physical Chemistry B 120.8 (2015): 1843–1848. Web.
3, 4, 6, 7, 12 | Conventional high-pressure ices III, IV,  VI, VII, and XII.
5      | Monoclinic ice V (testing).
16     | Negative-pressure ice XVI(16).  |Falenty, A., Hansen, T. C. & Kuhs, W. F. Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature 516, 231-233 (2014).
17     | Negative-pressure ice XVII(17).  |del Rosso, Leonardo, Milva Celli, and Lorenzo Ulivi. “Ice XVII as a Novel Material for Hydrogen Storage.” Challenges 8.1 (2017): 3.
0      | Hypothetical ice "0".  |Russo, J., Romano, F. & Tanaka, H. New metastable form of ice and its role in the homogeneous crystallization of water. Nat Mater 13, 733-739 (2014).
i      | Hypothetical ice "i".  |Fennell, C. J. & Gezelter, J. D. Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice I h. J. Chem. Theory Comput. 1, 662-667 (2005).
C0-II  | Filled ice C0 (Alias of 17). |Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560-3564 (2013).
C1     | Filled ice C1 (Alias of 2d).
C2     | Filled ice C2 (Alias of 1c).
sTprime | Filled ice sT' |Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560-3564 (2013).
CS1, CS2, TS1, HS1 | Clathrate hydrates CS1 (sI), CS2 (sII), TS1 (sIII), and HS1 (sIV).  |Matsumoto, M. & Tanaka, H. On the structure selectivity of clathrate hydrates. J. Phys. Chem. B 115, 8257-8265 (2011).
RHO    | Hypothetical ice at negative pressure ice 'sIII'. |Huang, Y et al. “A New Phase Diagram of Water Under Negative Pressure: the Rise of the Lowest-Density Clathrate S-III.” Science Advances 2.2 (2016): e1501010–e1501010.

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
