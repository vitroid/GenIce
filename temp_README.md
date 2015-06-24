# GenIce
Generate proton-disordered ice structures for GROMACS.

##Requirements
* Python 3
* NetworkX
* numpy

##Usage
%%usage%%
##Example
* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format:

    ./genice --cages --g12 co2 --g14 co2 --water tip4p CS1 > cs1.gro

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
