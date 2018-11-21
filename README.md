![Logo](https://raw.githubusercontent.com/vitroid/GenIce/develop/logo/genice-v0.png)
# GenIce

A Swiss army knife to generate hydrogen-disordered ice structures.

## Requirements

* Python 3
* NetworkX
* numpy
* svgwrite
* cif2ice
* countrings
* VPython (optional)

## Installation
GenIce is registered to [PyPI (Python Package Index)](https://pypi.python.org/pypi/GenIce). 
Install with pip3.

    pip3 install genice

## Uninstallation

    pip3 uninstall genice
    
## Usage

    usage: genice [-h] [--version] [--rep REP REP REP] [--dens DENS] [--seed SEED]
                  [--format gmeqdypoc] [--water model] [--guest D=empty]
                  [--Guest 13=me] [--Group 13=bu-:0] [--anion 3=Cl]
                  [--cation 3=Na] [--nodep] [--asis] [--debug] [--quiet]
                  Type
    
    GenIce is a swiss army knife to generate hydrogen-disordered ice structures.
    (version 0.22.11)
    
    positional arguments:
      Type                  Crystal type (1c,1h,etc. See
                            https://github.com/vitroid/GenIce for available ice
                            structures.)
    
    optional arguments:
      -h, --help            show this help message and exit
      --version, -V         show program's version number and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell in x,y, and z directions. [1,1,1]
      --dens DENS, -d DENS  Specify the ice density in g/cm3
      --seed SEED, -s SEED  Random seed [1000]
      --format gmeqdypoc, -f gmeqdypoc
                            Specify file format [g(romacs)|m(dview)|e(uler)|q(uate
                            rnion)|d(igraph)|y(aplot)|p(ython
                            module)|o(penScad)|c(entersofmass)|r(elative com)]
                            [gromacs]
      --water model, -w model
                            Specify water model. (tip3p, tip4p, etc.) [tip3p]
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
      --asis                Assumes all given HB pairs to be fixed. No shuffle and
                            no depolarization.
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.


Use `./genice.x` instead of `genice` if you want to use GenIce without installation. 

## Examples

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

        genice --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing
THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

        genice -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro

* With the aid of VPython, you can render and handle the ice structure
  directly in the web browser. (You must install VPython separately.)

        genice CS1 --format v

* You can read a .gro file as a unit cell of ice and convert to other format.  "Ow" and "Hw" are the atom name of water in the file.

        genice gromacs[mylattice.gro:Ow:Hw] --format scad --water tip5p > mylattice.scad

    (If you are using csh and descendants, you may need to escape `[` and `]` characters.)

## Basics

The program generates various ice lattice with proton disorder and without defect.  Total dipole moment is always set to zero (except the case you specify `--nodep` option).  The minimal structure (with --rep 1 1 1 option) is not always the unit cell of the lattice because it is difficult to deal with the hydrogen bond network topology of tiny lattice under periodic boundary condition.  Note that the generated structure is not optimal according to the potential energy.

* To get a large repetition of ice Ih in XYZ format,

        genice --rep 8 8 8 1h --format xyz > 1hx888.xyz

* To get a ice V lattice of different hydrogen order in CIF format, use `-s` option to specify the random seed.

        genice 5 -s 1024 --format cif > 5-1024.cif

* To obtain a ice VI lattice with different density and with TIP4P water model in gromacs format, use `--dens x` option to specify the density in g cm<sup>-3</sup>.

        genice 6 --dens 1.00 --format g --water tip4p > 6d1.00.gro

## File conversion

GenIce is a modular program; it reads a unit cell data from a lattice plugin defined in the lattices folder, put water and guest molecules using a molecule plugin defined in the molecules/ folder, and output in various formats using a format plugin defined in the formats/ folder. You can write your own plugins to extend GenIce. Some plugins also accept options. You can not only use GenIce to make a known ice structure but also use it to convert some file format to another.

* If you want to load a .gro file named "cs1.gro' and output it in yaplot format with using tip5p water model, just type the following.

        genice gromacs[cs1.gro:O:H] --format y --water tip5p --asis > cs1.yap

    where O and H are the atom names of water defined in the input .gro
file. You can use regular expression for hydrogen atom name. `--asis`
option avoids the network rearrangements. 

* If you want to let genice ignore hydrogen bonds and assume them from positions of oxygen atoms, specify the atom name of oxygen only.

        genice gromacs[cs1.gro:O] --format y > cs1.yap

* Some zeolites share the network topology with low-density ices. If you want to retrieve a zeolite ITT structure from [IZA structure database](http://www.iza-structure.org/databases) to prepare a low-density ice, try the following command:

        genice zeolite[ITT] -r 1 1 1 > ITT.gro



## Clathrate hydrates

For clathrate hydrates, you can prepare the lattice with cages partially occupied by various guest molecules.

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS
.gro format: (60% of small cages are filled with co2 and 40% are methane)

        genice -g 12=co2*0.6+me*0.4 -g 14=co2 --water tip4p CS1 > cs1.gro

* To make a CS2 clathrate hydrate structure of TIP5P water containing THF molecules in the large cage, while only one cage is filled with methane molecule, first just run genice without guest specifications:

        genice CS2 > CS2.gro
        
    The list of cages will be output as follows:

        INFO   Cage types: ['12', '16']
        INFO   Cage type 12: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183}
        INFO   Cage type 16: {136, 137, 138, 139, 140, 141, 142, 143, 16, 17, 18, 19, 20, 21, 22, 23, 160, 161, 162, 163, 164, 165, 166, 167, 40, 41, 42, 43, 44, 45, 46, 47, 184, 185, 186, 187, 188, 189, 190, 191, 64, 65, 66, 67, 68, 69, 70, 71, 88, 89, 90, 91, 92, 93, 94, 95, 112, 113, 114, 115, 116, 117, 118, 119}

    This indicates that there are two types of cages named `12` and `16`.  Fill the `16` cages with THF and put a methane molecule in the `0`th cage of type `12` as follows:
    
        genice CS2 -g 16=uathf -G 0=me > CS2.gro

Although only a few kinds of guest molecules are preset, you can easily prepare new guest molecules as a module. Here is an example for the ethlene oxide molecule.

`eo.py`

    import numpy as np
    # United-atom EO model with a dummy site
    LOC = 0.1436 # nm
    LCC = 0.1472 # nm
    
    Y = (LOC**2 - (LCC/2)**2)**0.5
    
    sites = np.array([[ 0.,    0., 0. ],
                      [-LCC/2, Y,  0. ],
                      [+LCC/2, Y,  0. ],])

    mass = np.array([16,14,14])
    # center of mass
    CoM = np.dot(mass, sites) / np.sum(mass)
    sites -= CoM
    
    atoms = ["O","C","C"]
    labels = ["Oe","Ce","Ce"]
    name = "EO"

Write the code in eo.py and put it in the user module folder at one of the following paths.

| OS | folder path |
|-------|-----------------------------|
|  | `./molecules`  |
| Linux | `~/.github/GenIce/molecules` |
| MacOS | `~/Library/Application Support/GenIce/molecules` |

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

    INFO   Hints:
    INFO     Cage types: ['12', '14', '15']
    INFO     Cage type 12: {0, 1, 2, 3, 4, 5, 14, 15, 16, 17, 18, 19, 28, 29, 30, 31, 32, 33, 42, 43, 44, 45, 46, 47, 56, 57, 58, 59, 60, 61, 70, 71, 72, 73, 74, 75, 84, 85, 86, 87, 88, 89, 98, 99, 100, 101, 102, 103}
    INFO     Cage type 14: {6, 7, 8, 9, 20, 21, 22, 23, 34, 35, 36, 37, 48, 49, 50, 51, 62, 63, 64, 65, 76, 77, 78, 79, 90, 91, 92, 93, 104, 105, 106, 107}
    INFO     Cage type 15: {10, 11, 12, 13, 24, 25, 26, 27, 38, 39, 40, 41, 52, 53, 54, 55, 66, 67, 68, 69, 80, 81, 82, 83, 94, 95, 96, 97, 108, 109, 110, 111}
    INFO     Cages adjacent to dopant 2: {9, 2, 28, 97}
    INFO     Cages adjacent to dopant 0: {9, 2, 28, 7}

It indicates that the nitrogen is surrounded by cages with ids 9, 2, 28, and 7.  Types for these cages can also be found in the info.  Then, we put the Bu- group (minus does not mean ions) in these cages adjacent dopant 0.

    genice HS1 -c 0=N -a 2=Br -H 9=Bu-:0 -H 2=Bu-:0 -H 28=Bu-:0 -H 7=Bu-:0 --nodep > HS1.gro

Here the option `-H` specifies the group by `-H (cage id)=(group name):(root)`, and root is the nitrogen that is specified by `-c` (cation) option.
 
 
### Placement of TBAB in the lattice module

*Under preparation*

It is more convenient if the lattice of the semiclathrate hydrate contains molecular ions in the appropriate locations in advance.  Here we explain the way to make the special module for semclathrates.

## AnalIce command

AnalIce is a variant of GenIce. AnalIce reads a Gromacs file and do not modify the molecular orientation or the hydrogen bond network topology.  AnalIce is prepared to use it for structure analysis.

For example, if you want to see the ring statistic of a given `.gro` file, use like this:


    analice input.gro -f _ring
	
If you want to replace water model from the original three-site one (described as OW, HW1, and HW2) to TIP4P-like four-site model, try

    analice input.gro -O OW -H HW[12] -w tip4p 

All the output formats are also available for AnalIce.


## Output formats

Name |Application | extension | water | solute | HB | remarks
-------|------------|-----------|----------|---------|-----|---
`cif`    |CIF         | `.cif`      | Atomic positions | Atomic positions | none |Experimental
`g`, `gromacs`      |[Gromacs](http://www.gromacs.org)     | `.gro`      | Atomic positions | Atomic positions | none|
`m`, `mdview`      |MDView      | `.mdv`      | Atomic positions | Atomic positions | auto|
`o`, `openscad`      |[OpenSCAD](http://www.openscad.org)    | `.scad`     | Center of mass | none | o | See tests/art/openscad for usage.
`povray`      |Povray | `.pov`     | Atomic positions | Atomic Positions | o | 
`towhee`      |TowHee    | `.coords`(?)      | Atomic positions | Atomic positions | none|
`xyz`    |XYZ         | `.xyz`      | Atomic positions | Atomic positions | none |Experimental
`exyz`    |[extended XYZ](http://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html)         | `.xyz`      | Atomic positions | Atomic positions | none |Extended XYZ format defined in Open Babel|
`exyz2`    |[extended XYZ](http://libatoms.github.io/QUIP/io.html#extendedxyz)         | `.xyz`      | Atomic positions | Atomic positions | none |Extended XYZ format defined in QUIP|
`y`, `yaplot`      |[Yaplot](https://github.com/vitroid/Yaplot)      | `.yap`      | Atomic positions | Atomic positions |o | It renders (1) HB paths to reduce the net polarization. (2) Rings in the structure. (3) Molecular configurations and the HB network.
`e`, `euler`      |Euler angles| `.nx3a`     | Rigid rotor | none | none|
`q`, `quaternion`      |Quaternions | `.nx4a`     | Rigid rotor | none |none|
`d`, `digraph`      |Digraph     | `.ngph`     | none | none | o |
`graph`  |Graph       | `.ngph`     | none | none | o | Experimental.
`c`, `com`      |CenterOfMass| `.ar3a`     | Center of mass | none | none |
`r`, `rcom`      |Relative CoM| `.ar3r`     | Center of mass | none | none | In fractional coordinate system.
`p`, `python`      |Python module | `.py`     | Center of mass | none | none | Under development.
`v`, `vpython`      |Direct visualization |     | Atomic positions | none | o | Display the structure in the browser using VPython.  You must install VPython separately. 
`svg_poly`      |SVG polygon| `.svg`     | Center of mass | none | o | See tests/art/svg for usage.
`_ring`      |Ring phase statistics |     | |  | | Statistical test suite 1: Check the appearance frequencies of the ring phases as a test for the intermediate-range disorder.
`_KG`      |Kirkwood G(r)|     | |  | | Statistical test suite 2: Calculate G(r) for checking long-range disorder in molecular orientations.

You can prepare your own file formats.  Write the module and put it in
one of the following paths.

| OS | folder path |
|-------|-----------------------------|
|  | `./formats`  |
| Linux | `~/.github/GenIce/formats` |
| MacOS | `~/Library/Application Support/GenIce/formats` |

Internally, there are seven stages to generate an ice structure.

1. Cell repetition.
2. Random graph generation and replication.
3. Apply ice rule.
4. Depolarize.
5. Determine orientations of the water molecules.
6. Place atoms in water molecules.
7. Place atoms in guests.

In the format plugin, you define the hook functions that are invoked after processing each stage. 

## Ice structures
<!-- references removed. -->

Symbol | Description
-------|------------
1h, 1c | Most popular Ice I (hexagonal or cubic)|
2|Proton-ordered ice II |
2d     | Hypothetical Proton-disordered Ice II.[Nakamura 2015]
3, 4, 6, 7, 12 | Conventional high-pressure ices III, IV,  VI, VII, and XII.[Lobban 1998]
5      | Monoclinic ice V (testing).
16     | Negative-pressure ice XVI(16).[Falenty 2014]
17     | Negative-pressure ice XVII(17).[del Rosso 2016]
0      | Hypothetical ice "0".[Russo 2014]
i      | Hypothetical ice "i". = Zeolite BCT.[Fennell 2005]
iceA, iceB   | Hypothetical ices "A" and "B".[Baez 1995]
quartzlike | Hypothetical "quartzlike" ice. [Svishchev 1996]
C0, C0-II  | Filled ice C0 (Alias of 17).[Smirnov 2013]
C1     | Filled ice C1 (Alias of 2).[Londono 1988]
C2     | Filled ice C2 (Alias of 1c).[Vos 1993]
sTprime | Filled ice "sT'" [Smirnov 2013]
CS1, CS2, CS4, TS1, HS1, HS2, HS3| Clathrate hydrates, Kosyakov's nomenclature. [Kosyakov 1999] 
sI, sII, sIII, sIV, sV, sVII, sH | Clathrate hydrates, Jeffrey's nomenclature. Jeffrey 1984]
RHO    | Hypothetical ice at negative pressure ice "sIII".[Huang 2016]
FAU    | Hypothetical ice at negative pressure ice "sIV". [Huang 2017]
CRN1, CRN2, CRN3 | 4-coordinated continuous random network [Mousseau 2005]
Struct01 .. Struct84 | Space Fullerenes [Dutour Sikiric 2010]
A15, sigma, Hcomp, Z, mu, zra-d, 9layers, 6layers, C36, C15, C14, delta, psigma | Space Fullerenes, Aliases of the Struct?? series.  See the data source for their names. [Dutour Sikiric 2010]
T      | Space fullerene type T,[Dutour Sikiric 2010] II+IVa. [Karttunen 2011]
xFAU[2], xFAU[4], xFAU[16], ... | Aeroices, i.e. extended FAU.[Matsui 2017]
iceR   | Partial plastic ice R [Mochizuki 2014].
iceT   | Partial plastic ice T [Hirata 2017].
iceT2  | Partial plastic ice T2 [Yagasaki 2018].
prism[4], prism[5], prism[6], ... | Ice nanotubes. [Koga 2001].
gromacs[filename]| Read a .gro file as a unit lattice of an ice.  See the output of `genice gromacs` for usage.  Note that only water molecules will be obtained. 
zeolite[XYZ]|Retrieve cif file of Zeolite XYZ from [IZA structure database](http://www.iza-structure.org/databases) as a unit lattice of an ice. Install [cif2ice](https://github.com/vitroid/cif2ice) separately to use it. (Experimental)
cif[filename]|Retrieve cif file as a unit lattice of an ice. Install [cif2ice](https://github.com/vitroid/cif2ice) separately to use it. (Experimental)

Ice names with double quotations are not experimentally verified.

You can prepare your own ice structures.  Write the module and put it in
one of the following paths.

| OS | folder path |
|-------|-----------------------------|
|  | `./lattices`  |
| Linux | `~/.github/GenIce/lattices` |
| MacOS | `~/Library/Application Support/GenIce/lattices` |

[cif2ice](https://github.com/vitroid/cif2ice) is a tool to retrieve a
cif file of zeolite from [IZA structure database](http://www.iza-structure.org/databases) and prepare a lattice
module in the path above.

Note: Some structures in different frameworks are identical.

CH/FI|CH  |ice|FK |Zeolite/Silica
-----|----|---|---|-----------
sI   |CS1 |-  |A15|MEP
sII  |CS2 |16 |C15|MTN
sIII |TS1 |-  |sigma|-  
sIV  |HS1 |-  |Z  |-  
sV   |HS2 |-  |*  |-  
sVII |CS4 |-  |*  |SOD
sH   |HS3 |-  |*  |DOH
C0   |-   |17 |*  |-  
C1   |-   |2  |*  |-  
C2   |-   |1c |*  |Cristobalite
     |-   |1h |*  |Tridymite
	 |-   |3  |*  |Keatite

FI: Filled ices; CH: Clathrate hydrates; FK:Frank-Kasper duals; Zeo: Zeolites.

-: No correspondence; *: Non-FK types.

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.
## Water models
A water model can be chosen with `--water` option.

symbol   | type
---------|-----
`3site`, `tip3p`  | 3-site TIP3P (default)
`4site`, `tip4p`  | 4-site TIP4P
`ice`             | TIP4P/ice
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


You can prepare your own guest molecules.  Write the module and put it in
one of the following paths.

| OS | folder path |
|-------|-----------------------------|
|  | `./molecules`  |
| Linux | `~/.github/GenIce/molecules` |
| MacOS | `~/Library/Application Support/GenIce/molecules` |

# References

* L. A. BÁEZ, P. CLANCY, J. Chem. Phys. 103, 9744–9755 (1998). [DOI: 10.1063/1.469938](http://doi.org/10.1063/1.469938)
* L. del Rosso, M. Celli, L. Ulivi, Nat Commun 2016, 7, 13394. [DOI: 10.1038/ncomms13394](http://dx.doi.org/10.1038/ncomms13394)
* M. Dutour Sikirić, O. Delgado-Friedrichs, M. Deza, Acta
  Crystallogr. A 2010, 66, 602. [DOI: 10.1107/S0108767310022932](http://dx.doi.org/10.1107/S0108767310022932)
* A. Falenty, T. C. Hansen, W. F. Kuhs, Nature 2014, 516, 231. [DOI: 10.1038/nature14014](http://dx.doi.org/10.1038/nature14014)
* C. J. Fennell, J. D. Gezelter, J. Chem. Theory Comput. 2005,
  1, 662. [DOI: 10.1021/ct050005s](http://dx.doi.org/10.1021/ct050005s)
* M. Hirata, T. Yagasaki, M. Matsumoto, H. Tanaka, to be published in
  Langmuir (2017). [DOI: 10.1021/acs.langmuir.7b01764](http://dx.doi.org/10.1021/acs.langmuir.7b01764)
* Y. Huang, C. Zhu, L. Wang, X. Cao, Y. Su, X. Jiang, S. Meng,
  J. Zhao, X. C. Zeng, Science Advances 2016, 2, e1501010. [DOI: 10.1126/sciadv.1501010](http://dx.doi.org/10.1126/sciadv.1501010)
* Y. Huang, C. Zhu, L. Wang, J. Zhao, X. C. Zeng,
  Chem. Phys. Lett. 2017, 671, 186. [DOI: 10.1016/j.cplett.2017.01.035](http://dx.doi.org/10.1016/j.cplett.2017.01.035)
* G. A. Jeffrey, In Inclusion Compounds; J. L. Atwood,
J. E. D. Davies, D. D. MacNicol, Eds.; Academic Press: London, 1984, Vol. 1, Chap. 5.
* A. J. Karttunen, T. F. Fässler, M. Linnolahti, T. A. Pakkanen, Inorg
  Chem 2011, 50, 1733. [DOI: 10.1021/ic102178d](http://dx.doi.org/10.1021/ic102178d)
* K Koga, GT Gao, H Tanaka, XC Zeng, Nature 412 (6849), 802-805 [DOI:10.1038/35090532](http://dx.doi.org/10.1038/35090532)
* V. I. Kosyakov, T. M. Polyanskaya, J. Struct. Chem. 1999, 40, 239. [DOI:10.1007/BF02903652](http://dx.doi.org/10.1007/BF02903652)
* C. Lobban, J. L. Finney, W. F. Kuhs, Nature 1998, 391, 268. [DOI: 10.1038/34622](http://dx.doi.org/10.1038/34622)
* D. Londono, W. F. Kuhs, J. L. Finney, Nature 1988, 332, 141. [DOI: 10.1038/332141a0](http://dx.doi.org/10.1038/332141a0)
* T. Matsui, M. Hirata, T. Yagasaki, M. Matsumoto, H. Tanaka,
  J. Chem. Phys. 147, 091101 (2017). [DOI: /10.1063/1.4994757](http://dx.doi.org/10.1063/1.4994757)
* K. Mochizuki, K. Himoto, and M. Matsumoto,
  Phys. Chem. Chem. Phys. 16, 16419–16425 (2014). [DOI: 10.1039/c4cp01616e](http://dx.doi.org/10.1039/c4cp01616e)
* N. Mousseau, G. T. Barkema, Curr. Opin. Solid State
  Mater. Sci. 2001, 5, 497. [DOI: 10.1016/S1359-0286(02)00005-0](https://doi.org/10.1016/S1359-0286%2802%2900005-0)
* T. Nakamura, M. Matsumoto, T. Yagasaki, H. Tanaka, J. Phys. Chem. B
  2015, 120, 1843. [DOI: 10.1021/acs.jpcb.5b09544](http://dx.doi.org/10.1021/acs.jpcb.5b09544)
* J. Russo, F. Romano, H. Tanaka, Nat Mater 2014, 13, 733. [DOI: 10.1038/NMAT3977](http://dx.doi.org/10.1038/NMAT3977)
* G. S. Smirnov, V. V. Stegailov, J Phys Chem Lett 2013, 4, 3560. [DOI: 10.1021/jz401669d](http://dx.doi.org/10.1021/jz401669d)
* I. M. Svishchev and P. G. Kusalik,Phys. Rev. B 53, R8815–R8817 (1996). [DOI: 10.1103/PhysRevB.53.R8815](http://doi.org/10.1103/PhysRevB.53.R8815)
* W. L. Vos, L. W. Finger, R. J. Hemley, H. Mao,
  Phys. Rev. Lett. 1993, 71, 3150. [DOI: 10.1103/PhysRevLett.71.3150](http://dx.doi.org/10.1103/PhysRevLett.71.3150)
* T. Yagasaki, M. Matsumoto, and H. Tanaka, J. Phys. Chem. B 122, 7718–7725 (2018). [DOI: 10.1021/acs.jpcb.8b04441](http://doi.org/10.1021/acs.jpcb.8b04441)

# The algorithm and how to cite it.

The algorithm to make a depolarized hydrogen-disordered ice is explained in our recent paper:

M. Masakazu, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered
Ice Generator",  J. Comput. Chem. 39, 61-64 (2017). [DOI: 10.1002/jcc.25077](http://dx.doi.org/10.1002/jcc.25077)

# How to contribute

GenIce is served at the GitHub (https://github.com/vitroid/GenIce) as an open source software since 2015. Feedbacks, proposals for improvements and extensions, and bug fixes are sincerely appreciated. Developers and test users are also welcome.
