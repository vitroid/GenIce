# GenIce
Generate proton-disordered ice structures for GROMACS.

##Requirements
* Python 3
* NetworkX
* numpy

##Usage
    usage: genice [-h] [--rep REP REP REP] [--dens DENS] [--seed SEED]
                  [--format gmeqd] [--water model] [--cages]
                  Type

    positional arguments:
      Type                  Crystal type (1c,1h,etc.)

    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [2,2,2]
      --dens DENS, -d DENS  Specify the ice density in g/cm3
      --seed SEED, -s SEED  Random seed [1000]
      --format gmeqd, -f gmeqd
                            Specify file format
                            [g(romacs)|m(dview)|e(uler)|q(uaternion)|d(igraph)]
      --water model, -w model
                            Specify water model.
      --cages, -c           Also output the cage positions. (g and m format only)

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
* Conventional ices 3, 4,  6, 7, and 12.
* Monoclinic ice 5 (testing).
* Ice 16.  [Falenty, A., Hansen, T. C. & Kuhs, W. F. Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature 516, 231–233 (2014).]
* Hypothetical ice 0.  [Russo, J., Romano, F. & Tanaka, H. New metastable form of ice and its role in the homogeneous crystallization of water. Nat Mater 13, 733–739 (2014).]
* Hypothetical ice i.  [Fennell, C. J. & Gezelter, J. D. Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice I h. J. Chem. Theory Comput. 1, 662–667 (2005).]
* Hypothetical ices sT' and C0-II.  [Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560–3564 (2013).]
* Clathrate hydrates CS1, CS2, HS1, and TS1.  [Matsumoto, M. & Tanaka, H. On the structure selectivity of clathrate hydrates. J. Phys. Chem. B 115, 8257–8265 (2011).]

Please ask vitroid@gmail.com to add new ice structures.
##Water models
* 3-site: TIP3P (default)
* 4-site: TIP4P
* 5-site: TIP5P
