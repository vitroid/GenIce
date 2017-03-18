# GenIce
Generate proton-disordered ice structures for GROMACS.

## Requirements

* Python 3
* NetworkX
* numpy

Note: WinPython includes all of these requirements.
## Installation
GenIce is registered to PyPI (Python Package Index). 
Install with pip3.

    pip3 install genice

## Uninstallation

    pip3 uninstall genice
    
## Usage
%%usage%%

## Example

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

        genice --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format:

        genice -g 12=co2 -g 14=co2 --water tip4p CS1 > cs1.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing
THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

        genice -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro

## Structure generation
The program generates various ice lattice with proton disorder and
without defect.  Total dipole moment is always set to zero.  The
minimal structure (with --rep 1 1 1 option) is not always the unit
cell of the lattice because it is difficult to deal with the hydrogen
bond network topology of tiny lattice under periodic boundary
condition.  Note that the generated structure is not optimal according
to the potential energy.
## Output formats

Format |Application | extension | water | solute | HB | remarks
-------|------------|-----------|----------|---------|-----|---
g      |[Gromacs](http://www.gromacs.org)     | .gro      | Atomic positions | Atomic positions | none|
m      |MDView      | .mdv      | Atomic positions | Atomic positions | auto|
e      |Euler angles| .nx3a     | Rigid rotor | none | none|
q      |Quaternions | .nx4a     | Rigid rotor | none |none|
d      |Digraph     | .ngph     | none | none | o |
y      |[Yaplot](https://github.com/vitroid/Yaplot)      | .yap      | o | o |none |
o      |[OpenSCAD](http://www.openscad.org)    | .scad     | Center of mass | none | o |
c      |CenterOfMass| .ar3a     | Center of mass | none | none |
p      |Python module | .py     | Center of mass | none | none | Under development.


## Ice structures

Symbol | Description| Remarks and data sources
-------|------------|----------
1h, 1c | Most popular Ice I (hexagonal or cubic)|
2d     | Hypothetical Proton-disordered Ice II. |Nakamura, Tatsuya et al. “Thermodynamic Stability of Ice II and Its Hydrogen-Disordered Counterpart: Role of Zero-Point Energy.” The Journal of Physical Chemistry B 120.8 (2015): 1843–1848. Web.
3, 4, 6, 7, 12 | Conventional high-pressure ices III, IV,  VI, VII, and XII.|
5      | Monoclinic ice V (testing).|
16     | Negative-pressure ice XVI(16).  |Falenty, A., Hansen, T. C. & Kuhs, W. F. Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature 516, 231-233 (2014).
17     | Negative-pressure ice XVII(17).  |del Rosso, Leonardo, Milva Celli, and Lorenzo Ulivi. “Ice XVII as a Novel Material for Hydrogen Storage.” Challenges 8.1 (2017): 3.
0      | Hypothetical ice "0".  |Russo, J., Romano, F. & Tanaka, H. New metastable form of ice and its role in the homogeneous crystallization of water. Nat Mater 13, 733-739 (2014).
i      | Hypothetical ice "i". = Zeolite BCT? |Fennell, C. J. & Gezelter, J. D. Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice I h. J. Chem. Theory Comput. 1, 662-667 (2005).
C0-II  | Filled ice C0 (Alias of 17). |Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560-3564 (2013).
C1     | Filled ice C1 (Alias of 2d).
C2     | Filled ice C2 (Alias of 1c).
sTprime | Filled ice "sT'" |Smirnov, G. S. & Stegailov, V. V. Toward Determination of the New Hydrogen Hydrate Clathrate Structures. J Phys Chem Lett 4, 3560-3564 (2013).
CS1, CS2, TS1, HS1 | Clathrate hydrates CS1 (sI), CS2 (sII), TS1 (sIII), and HS1 (sIV).  |Matsumoto, M. & Tanaka, H. On the structure selectivity of clathrate hydrates. J. Phys. Chem. B 115, 8257-8265 (2011).
RHO    | Hypothetical ice at negative pressure ice "sIII". |Huang, Y et al. “A New Phase Diagram of Water Under Negative Pressure: the Rise of the Lowest-Density Clathrate S-III.” Science Advances 2.2 (2016): e1501010–e1501010.
FAU    | Hypothetical ice at negative pressure ice "sIV". | “Prediction of a New Ice Clathrate with Record Low Density: a Potential Candidate as Ice XIX in Guest-Free Form.” “Prediction of a New Ice Clathrate with Record Low Density: a Potential Candidate as Ice XIX in Guest-Free Form.” sciencedirect.com. N.p., n.d. Web. 21 Feb. 2017.
CRN1,CRN2, CRN3 | 4-coordinated continuous random network| A model for low density amorphous ice. Mousseau, N, and G T Barkema. “Fast Bond-Transposition Algorithms for Generating Covalent Amorphous Structures.” Current Opinion in Solid State and Materials … 5.6 (2001): 497–502. Web.
Struct01 .. Struct84 | Space Fullerenes | Frank-Kasper type clathrate structures.    Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.
A15, sigma, Hcomp, Z, mu, zra-d, 9layers, 6layers, C36, C15, C14, delta, psigma | Space Fullerenes | Aliases of the Struct?? series.  See the data source for their names.  Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Ice names with double quotations are not experimentally verified.

Note: Some structures are identical.

Nomenclature |      |     |     |     |     |      |   |   |References
-------------|------|-----|-----|-----|-----|------|---|---|---
Frank-Kasper dual | A15 | C15 | sigma | Z   | C14| *| * | *| *: Non-FK types.  Frank, F.C., and JS Kasper. “Complex Alloy Structures Regarded as Sphere Packings. II. Analysis and Classification of Representative Structures.” Acta Crystallographica 12.7 (1959): 483–499.
ice | - | 16 |- |- |- | - |- |-
Jeffrey | sI | sII | sIII | sIV |  sV | sVI@|sVII | sH* | Jeffrey, G A. “Hydrate Inclusion Compounds.” Inclusion Compounds 1 (1984): 135–190.  (*) sH was not named by Jeffrey.
Kosyakov| CS1 | CS2 | TS1 | HS1 |  HS2 | CS3@|CS4| HS3|Kosyakov, Viktor I, and T M Polyanskaya. “Using Structural Data for Estimating the Stability of Water Networks in Clathrate and Semiclathrate Hydrates.” Journal of Structural Chemistry 40.2 (1999): 239–245.
Zeolite | MEP | MTN |  -|- | - | - | SOD |DOH |[New Database of Zeolite Structures](http://www.iza-structure.org/databases/)
@: not included in GenIce.

### Common structures between pure ices and hydrates
Nomenclature |  |   |   |  | Remarks and References
----|----|----|----|----| ---
ice | 1c | 2  |16 | 17 |
filled ice | C2 | C1 | sII | C0 |

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.
## Water models
A water model can be chosen with `--water` option.

symbol   | type|Reference
---------|--------|-------
`3site`, `tip3p`  | 3-site TIP3P (default)
`4site`, `tip4p`  | 4-site TIP4P
`5site`, `tip5p`  | 5-site TIP5P
`6site`, `NvdE`   | 6-site NvdE  | Nada, H.; van der Eerden, J. An Intermolecular Potential Model for the Simulation of Ice and Water Near the Melting Point: a Six-Site Model of H 2 O. J. Chem. Phys. 2003, 118, 7401.

## Guest molecules

symbol | type 
-------|---------
`co2`    | CO<sub>2</sub>
`uathf`  | United atom 5-site THF  
`g12`,`g14`,`g15`,`g16` | A monatomic dummy site

