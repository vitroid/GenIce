# GenIce
A Swiss army knife to generate proton-disordered ice structures.

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

    usage: genice [-h] [--rep REP REP REP] [--dens DENS] [--seed SEED]
                  [--format gmeqdypoc] [--water model] [--guest D=empty]
                  [--Guest 13=me] [--Group 13=bu-:0] [--anion 3=Cl]
                  [--cation 3=Na] [--nodep] [--debug] [--quiet]
                  Type
    
    positional arguments:
      Type                  Crystal type (1c,1h,etc. See
                            https://github.com/vitroid/GenIce for available ice
                            structures.)
    
    optional arguments:
      -h, --help            show this help message and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [2,2,2]
      --dens DENS, -d DENS  Specify the ice density in g/cm3
      --seed SEED, -s SEED  Random seed [1000]
      --format gmeqdypoc, -f gmeqdypoc
                            Specify file format [g(romacs)|m(dview)|e(uler)|q(uate
                            rnion)|d(igraph)|y(aplot)|p(ython
                            module)|o(penScad)|c(entersofmass)|r(elative com)]
      --water model, -w model
                            Specify water model. (tip3p, tip4p, etc.)
      --guest D=empty, -g D=empty
                            Specify guest(s) in the cage type. (D=empty,
                            T=co2*0.5+me*0.3, etc.)
      --Guest 13=me, -G 13=me
                            Specify guest in the specific cage. (13=me, 32=co2,
                            etc.)
      --Group 13=bu-:0, -H 13=bu-:0
                            Specify the group. (-H 13=bu-:0, etc.)
      --anion 3=Cl, -a 3=Cl
                            Specify a monatomic anion that replaces a water
                            molecule. (3=Cl, 39=F, etc.)
      --cation 3=Na, -c 3=Na
                            Specify a monatomic cation that replaces a water
                            molecule. (3=Na, 39=NH4, etc.)
      --nodep               No depolarization.
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.


## Examples

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

        genice --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing
THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

        genice -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro

## Basics

The program generates various ice lattice with proton disorder and without defect.  Total dipole moment is always set to zero (except the case you specify `--nodep` option).  The minimal structure (with --rep 1 1 1 option) is not always the unit cell of the lattice because it is difficult to deal with the hydrogen bond network topology of tiny lattice under periodic boundary condition.  Note that the generated structure is not optimal according to the potential energy.

* To get a large repetition of ice Ih in XYZ format,

        genice --rep 8 8 8 1h --format xyz > 1hx888.xyz

* To get a ice V lattice of different hydrogen order in CIF format, use `-s` option to specify the random seed.

        genice 5 -s 1024 --format cif > 5-1024.cif

* To obtain a ice VI lattice with different density and with TIP4P water model in gromacs format, use `--dens x` option to specify the density in g cm<sup>-3</sup>.

        genice 6 --dens 1.00 --format g --water tip4p > 6d1.00.gro

## Clathrate hydrates

For clathrate hydrates, you can prepare the lattice with cages partially occupied by various guest molecules.

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format: (60% of small cages are filled with co2 and 40% are methane)

        genice -g 12=co2*0.6+me*0.4 -g 14=co2 --water tip4p CS1 > cs1.gro

* To make a CS2 clathrate hydrate structure of TIP5P water containing THF molecules in the large cage, while only one cage is filled with methane molecule, first just run genice without guest specifications:

        genice CS2 > CS2.gro
        
    The list of cages will be output as follows:

   ```
INFO   Cage types: ['12', '16']
INFO   Cage type 12: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183}
INFO   Cage type 16: {136, 137, 138, 139, 140, 141, 142, 143, 16, 17, 18, 19, 20, 21, 22, 23, 160, 161, 162, 163, 164, 165, 166, 167, 40, 41, 42, 43, 44, 45, 46, 47, 184, 185, 186, 187, 188, 189, 190, 191, 64, 65, 66, 67, 68, 69, 70, 71, 88, 89, 90, 91, 92, 93, 94, 95, 112, 113, 114, 115, 116, 117, 118, 119}
```
    This indicates that there are two types of cages named `12` and `16`.  Fill the `16` cages with THF and put a methane molecule in the `0`th cage of type `12` as follows:
    
        genice CS2 -g 16=uathf -G 0=me > CS2.gro

Although only a few kinds of guest molecules are preset, you can easily prepare new guest molecules as a module. Here is an example for the ethlene oxide molecule.

`eo.py`

```
here is the code.
```

Write the code in eo.py and put it in the user module folder at either

* `~/.github/GenIce/molecules` Linux
* `~/Library/Application Support/GenIce/molecules` MacOS

*Note*: multiple occupancy is not implemented. If it is required, make a module of a virtual molecule that contains multiple molecules.

## Doping ions

Small ions may replace the host molecules.  In that case, you can use `-a` and `-c` options to replace the specified water molecules with anions and cations.

The following example replaces the `0`th water molecule (in the replicated lattice) with Na cation and `1`st water molecule with Cl anion.  The hydrogen bonds around the ions are organized appropriately.

    genice CS2 --nodep -c 0=Na -a 1=Cl > CS2.gro

*Note 1*: The numbers of cations and anions must be the same.  Otherwise, ice rule is never satisfied and the program does not stop.  

*Note 2*: The option `--nodep` is also required because it is impossible to depolarize the structure containing ions.

*Note 3*: Protonic defects (H<sub>3</sub>O<sup>+</sup> and OH<sup>-</sup>) are not yet implemented.

## Semiclathrate hydrates

### Placement of a tetrabutylammonium ion

Let us assume that the id of the water molecule to be replaced by nitrogen of the TBA as zero.  Place the nitrogen as a cation and also replace the water 2 by the counterion Br.

    genice HS1 -c 0=N -a 2=Br --nodep > HS1.gro

Then you will see the following info.

```
INFO   Hints:
INFO     Cage types: ['12', '14', '15']
INFO     Cage type 12: {0, 1, 2, 3, 4, 5, 14, 15, 16, 17, 18, 19, 28, 29, 30, 31, 32, 33, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 70, 71, 72, 73, 74, 75, 84, 85, 86, 87, 88, 89, 98, 99, 100, 101, 102, 103}
INFO     Cage type 14: {6, 7, 8, 9, 20, 21, 22, 23, 34, 35, 36, 37, 48, 49, 50, 51, 62, 63, 64, 65, 76, 77, 78, 79, 90, 91, 92, 93, 104, 105, 106, 107}
INFO     Cage type 15: {10, 11, 12, 13, 24, 25, 26, 27, 38, 39, 40, 41, 52, 53, 54, 55, 66, 67, 68, 69, 80, 81, 82, 83, 94, 95, 96, 97, 108, 109, 110, 111}
INFO     Cages adjacent to dopant 2: {9, 2, 28, 97}
INFO     Cages adjacent to dopant 0: {9, 2, 28, 7}
```

It indicates that the nitrogen is surrounded by cages with ids 9, 2, 28, and 7.  Types for these cages can also be found in the info.  Then, we put the Bu- group (minus does not mean ions) in these cages adjacent dopant 0.

    genice HS1 -c 0=N -a 2=Br -H 9=Bu-:0 -H 2=Bu-:0 -H 28=Bu-:0 -H 7=Bu-:0 --nodep > HS1.gro

Here the option `-H` specifies the group by `-H (cage id)=(group name):(root)`, and root is the nitrogen that is specified by `-c` (cation) option.
 
 
### Placement of TBAB in the lattice module

*Under preparation*

It is more convenient if the lattice of the semiclathrate hydrate contains molecular ions in the appropriate locations in advance.  Here we explain the way to make the special module for semclathrates.

## Output formats

Name |Application | extension | water | solute | HB | remarks
-------|------------|-----------|----------|---------|-----|---
`g`, `gromacs`      |[Gromacs](http://www.gromacs.org)     | `.gro`      | Atomic positions | Atomic positions | none|
`m`, `mdview`      |MDView      | `.mdv`      | Atomic positions | Atomic positions | auto|
`e`, `euler`      |Euler angles| `.nx3a`     | Rigid rotor | none | none|
`q`, `quaternion`      |Quaternions | `.nx4a`     | Rigid rotor | none |none|
`d`, `digraph`      |Digraph     | `.ngph`     | none | none | o |
`graph`  |Graph       | `.ngph`     | none | none | o | Experimental.
`y`, `yaplot`      |[Yaplot](https://github.com/vitroid/Yaplot)      | `.yap`      | o | o |none |
`o`, `openscad`      |[OpenSCAD](http://www.openscad.org)    | `.scad`     | Center of mass | none | o |
`c`, `com`      |CenterOfMass| `.ar3a`     | Center of mass | none | none |
`r`, `rcom`      |Relative CoM| `.ar3r`     | Center of mass | none | none | In fractional coordinate system.
`p`, `python`      |Python module | `.py`     | Center of mass | none | none | Under development.
`cif`    |CIF         | `.cif`      | Atomic positions | Atomic positions | none |Experimental

## Ice structures
<!-- rreferences removed. -->

Symbol | Description
-------|------------
1h, 1c | Most popular Ice I (hexagonal or cubic)|
2|Proton-ordered ice II |
2d     | Hypothetical Proton-disordered Ice II.
3, 4, 6, 7, 12 | Conventional high-pressure ices III, IV,  VI, VII, and XII.
5      | Monoclinic ice V (testing).|
16     | Negative-pressure ice XVI(16).
17     | Negative-pressure ice XVII(17).
0      | Hypothetical ice "0".
i      | Hypothetical ice "i". = Zeolite BCT
C0, C0-II  | Filled ice C0 (Alias of 17).
C1     | Filled ice C1 (Alias of 2).
C2     | Filled ice C2 (Alias of 1c).
sTprime | Filled ice "sT'"
CS1, CS2, CS4, TS1, HS1, HS2, HS3, sI, sII, sIII, sIV, sV, sVII, sH | Clathrate hydrates
RHO    | Hypothetical ice at negative pressure ice "sIII".
FAU    | Hypothetical ice at negative pressure ice "sIV".
CRN1,CRN2, CRN3 | 4-coordinated continuous random network
Struct01 .. Struct84 | Space Fullerenes
A15, sigma, Hcomp, Z, mu, zra-d, 9layers, 6layers, C36, C15, C14, delta, psigma | Space Fullerenes, Aliases of the Struct?? series.  See the data source for their names.

Ice names with double quotations are not experimentally verified.

Note: Some structures in different frameworks are identical.

CH/FI|CH  |ice|FK |Zeo|
-----|----|---|---|---|
sI   |CS1 |-  |A15|MEP|
sII  |CS2 |16 |C15|MTN|
sIII |TS1 |-  |σ  |-  |
sIV  |HS1 |-  |Z  |-  |
sV   |HS2 |-  |*  |-  |
sVII |CS4 |-  |*  |SOD|
sH   |HS3 |-  |*  |DOH|
C0   |-   |17 |*  |-  |
C1   |-   |2  |*  |-  |
C2   |-   |1c |*  |-  |

FI: Filled ices; CH: Clathrate hydrates; FK:Frank-Kasper duals; Zeo: Zeolites.

-: No correspondence; *: Non-FK types; @: Not included in GenIce.

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.
## Water models
A water model can be chosen with `--water` option.

symbol   | type
---------|-----
`3site`, `tip3p`  | 3-site TIP3P (default)
`4site`, `tip4p`  | 4-site TIP4P
`5site`, `tip5p`  | 5-site TIP5P
`6site`, `NvdE`   | 6-site NvdE

## Guest molecules

symbol | type 
-------|---------
`co2`    | CO<sub>2</sub>
`me`     | United atom monatomic methane
`uathf`  | United atom 5-site THF  
`g12`,`g14`,`g15`,`g16` | A monatomic dummy site
`empty`  | Leave the cage empty.

