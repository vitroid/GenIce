![Logo](https://raw.githubusercontent.com/vitroid/GenIce/develop/logo/genice-v0.png)
# GenIce

A Swiss army knife to generate hydrogen-disordered ice structures.

version 1.0.7

## Requirements

* networkx>=2
* countrings>=0.1.7
* pairlist>=0.2.3
* yaplotlib>=0.1
* numpy

## Installation
GenIce is registered to [PyPI (Python Package Index)](https://pypi.python.org/pypi/GenIce). 
Install with pip3.

    pip3 install genice

## Uninstallation

    pip3 uninstall genice
    
## Usage

    usage: genice [-h] [--version] [--rep REP REP REP] [--dens DENS]
                  [--add_noise percent] [--seed SEED] [--format name]
                  [--water model] [--guest D=empty] [--Guest 13=me]
                  [--Group 13=bu-:0] [--anion 3=Cl] [--cation 3=Na]
                  [--visual visual] [--nodep] [--asis] [--debug] [--quiet]
                  Type
    
    GenIce is a swiss army knife to generate hydrogen-disordered ice structures.
    (version 1.0.7)
    
    positional arguments:
      Type                  Crystal type (1c, 1h, etc. See 
                            https://github.com/vitroid/GenIce for available ice 
                            structures.)
                            If you want to analyze your own structures, please try 
                            analice tool.
                             
                             
                            [Available lattice structures]
                             
                            1. Lattice structures served with GenIce
                             
                            0, ice0         Metastable ice "0".
                            11              Ice XI.
                            12, XII, ice12  Ice XII.
                            13, XIII, ice13 Ice XIII.
                            16, CS2, MTN, XVI, ice16, sII   Ice XVI.
                            17, C0-II, C0, XVII, ice17      Ice XVII.
                            1c, C2, Ic, ice1c               Ice Ic.
                            1h, Ih, ice1h   Ice Ih.
                            2, C1, II, ice2 Ice II; Hydrogen hydrate C1.
                            2D3             Trilayer honeycomb ice.
                            2d, ice2d       A hydrogen-disordered counterpart of 
                                            ice II.
                            3, III, ice3    Ice III.
                            4, IV, ice4     Ice IV.
                            4R              Orthogonalized ice IV.
                            5, V, ice5      Ice V.
                            5R              Orthogonalized ice V.
                            6, VI, ice6     Ice VI.
                            6h              Half lattice of ice VI.
                            7, VII, ice7    Ice VII.
                            8, VIII, ice8   Ice VIII.
                            9, IX, ice9     Ice IX.
                            A, iceA         Hypothetical ice A.
                            A15, Struct33   Cubic Structure I of clathrate hydrate.
                            B, iceB         Hypothetical ice B.
                            ice11_19        A candidate for an antiferroelectric 
                                            Ice XI #19.
                            ice2rect        Orthogonalized Ice II.
                            ----
                            (Undocumented) 1h_unit C14 C15 C36 CRN1 CRN2 CRN3 CS1 
                            CS4 DOH EMT FAU FK6layers FK9layers HS1 HS2 HS3 Hcomp 
                            Kcomp LTA MEP RHO SOD Struct01 Struct02 Struct03 
                            Struct04 Struct05 Struct06 Struct07 Struct08 Struct09 
                            Struct10 Struct11 Struct12 Struct13 Struct14 Struct15 
                            Struct16 Struct17 Struct18 Struct19 Struct20 Struct21 
                            Struct22 Struct23 Struct24 Struct25 Struct26 Struct27 
                            Struct28 Struct29 Struct30 Struct31 Struct32 Struct34 
                            Struct35 Struct36 Struct37 Struct38 Struct39 Struct40 
                            Struct41 Struct42 Struct43 Struct44 Struct45 Struct46 
                            Struct47 Struct48 Struct49 Struct50 Struct51 Struct52 
                            Struct53 Struct54 Struct55 Struct56 Struct57 Struct58 
                            Struct59 Struct60 Struct61 Struct62 Struct63 Struct64 
                            Struct65 Struct66 Struct67 Struct68 Struct69 Struct70 
                            Struct71 Struct72 Struct73 Struct74 Struct75 Struct76 
                            Struct77 Struct78 Struct79 Struct80 Struct81 Struct82 
                            Struct83 Struct84 T TS1 Z delta dtc i ice1h_unit iceR 
                            iceT iceT2 mu prism psigma sH sI sIII sIV sTprime sV 
                            sVII sigma xFAU xFAU2 zra-d
                             
                             
                            2. Lattice structures served by plugins
                             
                            (None)
                            ----
                             
                             
                            3. Lattice structures served locally
                             
                            (None)
                            ----
                             
                             
    
    optional arguments:
      -h, --help            show this help message and exit
      --version, -V         show program's version number and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell along a, b, and c axes. [1,1,1]
      --dens DENS, -d DENS  Specify the ice density in g/cm3 (Guests are not
                            included.)
      --add_noise percent   Add a Gauss noise with given width (SD) to the
                            molecular positions of water. The value 1 corresponds
                            to 1 percent of the molecular diameter of water.
      --seed SEED, -s SEED  Random seed [1000]
      --format name, -f name
                            Specify the output file format. [gromacs]
                             
                             
                            [Available formatters]
                             
                            1. Formatters served with GenIce
                             
                            _KG             Kirkworrd G factor.
                            _ringstat       Bond direction statistics.
                            d, digraph      Directed graph of HBs.
                            e, euler        Rigid rotor (Euler angle).
                            exmol           Extended XMol file format.
                            exyz            Extended XYZ format.
                            g, gromacs      Gromacs .gro file.
                            graph           Undirected graph of HBs.
                            m, mdview       MDView file (in Angdtrom).
                            mdv_au          MDView file (in au).
                            o, openscad     OpenSCAD.
                            p, python, reshape              Cell-reshaper.
                            povray          Povray.
                            q, quaternion   Rigid rotor (Quaternion).
                            rings           Show rings in Yaplot.
                            y, yaplot       Yaplot.
                            ----
                            (Undocumented) bdl c cif cif2 com r rcom towhee xyz
                             
                             
                            2. Formatters served by plugins
                             
                            _RDF            Radial Distribution Functions.
                            ----
                             
                             
                            3. Formatters served locally
                             
                            (None)
                            ----
                             
                             
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
      --visual visual       Specify the yaplot file to store the depolarization
                            paths. [""]
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


## Basics

The program generates various ice lattice with proton disorder and without defect.  Total dipole moment is always set to zero (except the case you specify `--nodep` option).  The minimal structure (with --rep 1 1 1 option) is not always the unit cell of the lattice because it is difficult to deal with the hydrogen bond network topology of tiny lattice under periodic boundary condition.  Note that the generated structure is not optimal according to the potential energy.

* To get a large repetition of ice Ih in XYZ format,

        genice --rep 8 8 8 1h --format xyz > 1hx888.xyz

* To get a ice V lattice of different hydrogen order in CIF format, use `-s` option to specify the random seed.

        genice 5 -s 1024 --format cif > 5-1024.cif

* To obtain a ice VI lattice with different density and with TIP4P water model in gromacs format, use `--dens x` option to specify the density in g cm<sup>-3</sup>.

        genice 6 --dens 1.00 --format g --water tip4p > 6d1.00.gro

GenIce is a modular program; it reads a unit cell data from a lattice plugin defined in the lattices folder, put water and guest molecules using a molecule plugin defined in the molecules/ folder, and output in various formats using a format plugin defined in the formats/ folder. You can write your own plugins to extend GenIce. Some plugins also accept options.

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

Write the code in eo.py. Make a folder named `molecules` in the current working directory and put it in.

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


    analice input.gro -f _ringstat
	
If you want to replace water model from the original three-site one (described as OW, HW1, and HW2) to TIP4P-like four-site model, try

    analice input.gro -O OW -H HW[12] -w tip4p 

All the output formats are also available for AnalIce.

### More examples

Load every 10 frames from a set of .gro files and output ring statistics in separate files.

    analice '%05d.gro' --framerange 0:1000000:10 -O OW -H HW[12] --format _ringstat -o '%04d.rstat'
	
Make V-structures (removal of quick librational motion of water) from the given set of .gro files.

    analice '%05d.gro' -O OW -H HW[12] -w tip3p --avgspan 25 > vstruct.gro

### Usage of `analice`

    usage: analice [options]
    
    GenIce is a swiss army knife to generate hydrogen-disordered ice structures.
    (version 1.0.7)
    
    positional arguments:
      File                  Input file(s). File type is estimated from the suffix. 
                            Files of different types cannot be read at a time. File
                             type can be specified explicitly with -s option.
                             
                             
                            [Available input file types]
                             
                            1. File types served with GenIce
                             
                            exyz            Extended XYZ format.
                            gro             Gromacs .gro file.
                            mdv             MDView file (in Angdtrom).
                            mdva            MDView file (in au).
                            ----
                            (Undocumented) ar3r nx3a
                             
                             
                            2. File types served by plugins
                             
                            (None)
                            ----
                             
                             
                            3. File types served locally
                             
                            (None)
                            ----
                             
                             
    
    optional arguments:
      -h, --help            show this help message and exit
      --version, -V         show program's version number and exit
      --format name, -f name
                            Specify the output file format. [gromacs]
                             
                             
                            [Available formatters]
                             
                            1. Formatters served with GenIce
                             
                            _KG             Kirkworrd G factor.
                            _ringstat       Bond direction statistics.
                            d, digraph      Directed graph of HBs.
                            e, euler        Rigid rotor (Euler angle).
                            exmol           Extended XMol file format.
                            exyz            Extended XYZ format.
                            g, gromacs      Gromacs .gro file.
                            graph           Undirected graph of HBs.
                            m, mdview       MDView file (in Angdtrom).
                            mdv_au          MDView file (in au).
                            o, openscad     OpenSCAD.
                            p, python, reshape              Cell-reshaper.
                            povray          Povray.
                            q, quaternion   Rigid rotor (Quaternion).
                            rings           Show rings in Yaplot.
                            y, yaplot       Yaplot.
                            ----
                            (Undocumented) bdl c cif cif2 com r rcom towhee xyz
                             
                             
                            2. Formatters served by plugins
                             
                            _RDF            Radial Distribution Functions.
                            ----
                             
                             
                            3. Formatters served locally
                             
                            (None)
                            ----
                             
                             
      --output %04d.gro, -o %04d.gro
                            Output in separate files.
      --water model, -w model
                            Replace water model. (tip3p, tip4p, etc.) [tip3p]
      --oxygen OW, -O OW    Specify atom name of oxygen in input Gromacs file.
                            ("O")
      --hydrogen HW[12], -H HW[12]
                            Specify atom name (regexp) of hydrogen in input
                            Gromacs file. ("H")
      --suffix gro, -s gro  Specify the file suffix explicitly. ((None)
      --filerange [from:]below[:interval]
                            Specify the number range for the input filename.
                            ("0:1000000")
      --framerange [from:]below[:interval]
                            Specify the number range for the input frames.
                            ("0:1000000")
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.
      --add_noise percent   Add a Gauss noise with given width (SD) to the
                            molecular positions of water. The value 1 corresponds
                            to 1 percent of the molecular diameter of water.
      --avgspan 1, -v 1     Output mean atomic positions of a given time span so
                            as to remove fast librational motions and to make a
                            smooth video. The values 0 and 1 specify no averaging.

## Output formats (`genice` and `analice`)

Name |Application | extension | water | solute | HB | remarks
-------|------------|-----------|----------|---------|-----|---
`cif, cif2 |CIF         | `.cif`      | Atomic positions | Atomic positions | none |Experimental
`g`, `gromacs`      |[Gromacs](http://www.gromacs.org)     | `.gro`      | Atomic positions | Atomic positions | none| Default format.
`m`, `mdview`      |MDView      | `.mdv`      | Atomic positions | Atomic positions | auto|
`mdv_au`      |MDView      | `.mdv`      | Atomic positions | Atomic positions | auto| In atomic unit.
`o`, `openscad`      |[OpenSCAD](http://www.openscad.org)    | `.scad`     | Center of mass | none | o | See tests/art/openscad for usage.
`povray`      |Povray | `.pov`     | Atomic positions | Atomic Positions | o | 
`towhee`      |TowHee    | `.coords`(?)      | Atomic positions | Atomic positions | none|
`xyz`    |XYZ         | `.xyz`      | Atomic positions | Atomic positions | none |Experimental
`exyz`    |[extended XYZ](http://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html)         | `.xyz`      | Atomic positions | Atomic positions | none |Extended XYZ format defined in Open Babel|
`exyz2`    |[extended XYZ](http://libatoms.github.io/QUIP/io.html#extendedxyz)         | `.xyz`      | Atomic positions | Atomic positions | none |Extended XYZ format defined in QUIP|
`y`, `yaplot`      |[Yaplot](https://github.com/vitroid/Yaplot)      | `.yap`      | Atomic positions | Atomic positions |o | It renders molecular configurations and the HB network.
`e`, `euler`      |Euler angles| `.nx3a`     | Rigid rotor | none | none|
`q`, `quaternion`      |Quaternions | `.nx4a`     | Rigid rotor | none |none|
`d`, `digraph`      |Digraph     | `.ngph`     | none | none | o |
`graph`  |Graph       | `.ngph`     | none | none | o | Experimental.
`c`, `com`      |CenterOfMass| `.ar3a`     | Center of mass | none | none |
`r`, `rcom`      |Relative CoM| `.ar3r`     | Center of mass | none | none | In fractional coordinate system.
`p`, `python`, `reshape`      |Python module | `.py`     | Center of mass | none | none | Under development.
`_ringstat`      |Ring phase statistics |     | |  | | Statistical test suite 1: Check the appearance frequencies of the ring phases as a test for the intermediate-range disorder.
`rings`      |[Yaplot](https://github.com/vitroid/Yaplot)      | `.yap`      | center of mass | none |o | It renders HB rings.
`_KG`      |Kirkwood G(r)|     | |  | | Statistical test suite 2: Calculate G(r) for checking long-range disorder in molecular orientations.

You can prepare your own file formats. Create a folder named `formats` in the current working directory and put the plugins in it. GenIce 1.0 no longer refers the files in `~/.genice` folder.

Internally, there are seven stages to generate an ice structure.

1. Cell repetition.
2. Random graph generation and replication.
3. Apply ice rule.
4. Depolarize.
5. Determine orientations of the water molecules.
6. Place atoms in water molecules.
7. Place atoms in guests.

In the format plugin, you define the hook functions that are invoked after processing each stage. 

## Ice structures (`genice` only)
<!-- references removed. -->

Symbol | Description
-------|------------
1h, ice1h, Ih   | Most popular Ice I (hexagonal)
1c, ice1c, Ic   | Cubic type of ice I
2, ice2, II     | Hydrogen-ordered ice II
2d, ice2d       | Hypothetical Hydrogen-disordered Ice II.[Nakamura 2015]
3, ice3, III    | Conventional high-pressure ice III.[Lobban 1998]
4, ice4, IV     | Metastable high-pressure ice IV.[Lobban 1998]
4R              | Ice IV with orthogonal unit cell. (testing)
5, ice5, V      | Monoclinic ice V (testing).
5R              | Ice V with orthogonal unit cell. (testing)
6, ice6, VI     | Conventional high-pressure ice VI.[Lobban 1998]
6h              | Half lattice of ice IV.
7, ice7, VII    | Conventional high-pressure ice VII.[Lobban 1998]
8, ice8, VIII   | Ice VIII, a hydrogen-ordered counterpart of ice VII.[Kuhs 1998]
9, ice9, IX     | Ice IX, a hydrogen-ordered counterpart of ice III.[Londono 1993]
ice11_19        | A candidate for an antiferroelectric ice 11; ice 11 type 19 in Ref. [Fan 2010]
12, ice12, XII  | Metastable high-pressure ice XII.[Lobban 1998]
13, ice13, XIII | Ice XIII, a hydrogen-ordered counterpart of ice V.[Salzmann 2006]
16, ice16, XVI  | Negative-pressure ice XVI.[Falenty 2014]
17, ice17, XVII | Negative-pressure ice XVII.[del Rosso 2016]
0, ice0         | Hypothetical ice "0".[Russo 2014]
i               | Hypothetical ice "i". = Zeolite BCT.[Fennell 2005]
A, iceA         | Hypothetical hydrogen-ordered ices "A" and "B".[Baez 1998]
B, iceB         | Hypothetical hydrogen-ordered ices "A" and "B".[Baez 1998]
C0, C0-II       | Filled ice C0 (Alias of 17).[Smirnov 2013]
C1              | Filled ice C1 (Alias of 2).[Londono 1988]
C2              | Filled ice C2 (Alias of 1c).[Vos 1993]
sTprime         | Filled ice "sT'". [Smirnov 2013]
CS1, CS2, CS4, TS1, HS1, HS2, HS3| Clathrate hydrates, Kosyakov's nomenclature. [Kosyakov 1999] 
sI, sII, sIII, sIV, sV, sVII, sH | Clathrate hydrates, Jeffrey's nomenclature. [Jeffrey 1984]
RHO             | Hypothetical ice at negative pressure ice "sIII".[Huang 2016]
FAU             | Hypothetical ice at negative pressure ice "sIV". [Huang 2017]
EMT             | Hypothetical ice with a large cavity.[Liu 2019]
DOH,MEP,MTN,SOD | Aliases of HS3, CS1, CS2, and CS4, respectively.
CRN1,CRN2,CRN3  | 4-coordinated continuous random network [Mousseau 2005]
Struct01 .. Struct84 | Space Fullerenes [Dutour Sikiric 2010]
A15, sigma, Hcomp, Kcomp, Z, mu, zra-d, FK9layers, FK6layers, C36, C15, C14, delta, psigma | Space Fullerenes, Aliases of the Struct?? series.  See the data source for their names. [Dutour Sikiric 2010]
T      | Space fullerene type T,[Dutour Sikiric 2010] II+IVa. [Karttunen 2011]
xFAU[2], xFAU[4], xFAU[16], ... | Aeroices, i.e. extended FAU.[Matsui 2017]
xFAU2[2], xFAU2[4], xFAU2[16], ... | Aeroices, i.e. extended FAU.[Matsui 2017] (Hydrogen bond orientations are modified.)
iceR   | Partial plastic ice R [Mochizuki 2014].
iceT   | Partial plastic ice T [Hirata 2017].
iceT2  | Partial plastic ice T2 [Yagasaki 2018].
dtc    | Ultralow-density ice containing cylindrical pores. [Matsui 2019]
prism[4], prism[5], prism[6], ... | Ice nanotubes. [Koga 2001].

Ice names with double quotations are not experimentally verified.

You can prepare your own ice structures. Create a folder named `lattices` in the current working directory and put the plugins in it. GenIce 1.0 no longer refers the files in `~/.genice` folder.

[cif2ice](https://github.com/vitroid/cif2ice) is a tool to retrieve a
cif file of zeolite from [IZA structure database](http://www.iza-structure.org/databases) and prepare a lattice
module in the path above.

Note: Some structures in different frameworks are identical.

CH/FI|CH  |ice|FK |Zeo|
-----|----|---|---|---|
sI   |CS1 |-  |A15|MEP|
sII  |CS2 |16 |C15|MTN|
sIII |TS1 |-  |sigma|-  |
sIV  |HS1 |-  |Z  |-  |
sV   |HS2 |-  |*  |-  |
sVII |CS4 |-  |*  |SOD|
sH   |HS3 |-  |*  |DOH|
C0   |-   |17 |*  |-  |
C1   |-   |2  |*  |-  |
C2   |-   |1c |*  |-  |

FI: Filled ices; CH: Clathrate hydrates; FK:Frank-Kasper duals; Zeo: Zeolites.

-: No correspondence; *: Non-FK types.

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.
## Water models (`genice` and `analice`)
A water model can be chosen with `--water` option.

symbol   | type
---------|-----
`3site`, `tip3p`  | 3-site TIP3P (default)
`4site`, `tip4p`  | 4-site TIP4P
`ice`             | TIP4P/ice
`5site`, `tip5p`  | 5-site TIP5P
`6site`, `NvdE`   | 6-site NvdE

## Guest molecules (`genice` only)

symbol | type 
-------|---------
`co2`    | CO<sub>2</sub>
`me`     | United atom monatomic methane
`uathf`  | United atom 5-site THF  
`g12`,`g14`,`g15`,`g16` | A monatomic dummy site
`empty`  | Leave the cage empty.


You can prepare your own guest molecules.  Create a folder named `molecules` in the current working directory and put the plugins in it. GenIce 1.0 no longer refers the files in `~/.genice` folder.

## Input files (`analice` only)

suffix | type
-------|-------------------
`gro`  | Gromacs
`mdv`  | mdview (Angstrom)
`mdva` | mdview (au)
`nx3a` | Rigid rotors (with euler angles)

You can prepare your own file loaders.  Create a folder named `loaders` in the current working directory and put the plugins in it. The plugin name is refered as the suffix of the file. (e.g. prepare pdb.py to load a *.pdb file.)

# Extra plugins
(New in v1.0)

Some extra plugins are available via python package index using pip command.

For example, you can install RDF plugin by the following command,

    % pip install genice-rdf
	
And use it as an output format to get the radial distribution functions.

    % genice TS1 -f _RDF > TS1.rdf.txt



## Output and analysis plugins
Analysis plugin is a kind of output plugin (specified with -f option). They are useful with analice command.

| pip name | GenIce option | Description | output format | requirements |
|----------|-------|-------------|---------------|--------------|
|[`genice-cage`](https://github.com/vitroid/genice-cage)|`-f _cage`| Detect cages and quasi-polyhedra (vitrites). | text, json | `countrings` |
|[`genice-rdf`](https://github.com/vitroid/genice-rdf)|`-f _RDF`| Radial distribution functions. | text |  |
|[`genice-svg`](https://github.com/vitroid/genice-svg)|`-f svg`<br />`-f png` | 2D graphics in SVG format.<br /> ... in PNG format.| SVG<br />PNG | `svgwrite` |
|[`genice-vpython`](https://github.com/vitroid/genice-vpython)|`-f vpython`| Display the structure in the browser using VPython.| (none) | `vpython` |
|[`genice-twist`](https://github.com/vitroid/genice-twist)|`-f twist`| Calculate the twist order parameter (and visualize) [Matsumoto 2019]| text (@BTWC)<br />SVG<br />PNG <br />yaplot | `twist-op`, `genice-svg` |

## Input plugins
Input plugins (a.k.a. lattice plugins) construct a crystal structure on demand.

| pip name   | GenIce usage    | Description  |requirements |
|------------|-----------------|--------------|-------------|
|[`genice-cif`](https://github.com/vitroid/genice-cif)| `genice cif[ITT.cif]`<br /> `genice zeolite[ITT]`| Read a local CIF file as an ice structure.<br />Read a structure from Zeolite DB. | `cif2ice` |

# References

* L.A. Báez, P. Clancy, 
  J. Chem. Phys. 103, 9744–9755 (1998).
  [DOI: 10.1063/1.469938](http://doi.org/10.1063/1.469938)
* L. del Rosso, M. Celli, L. Ulivi,
  Nat Commun 2016, 7, 13394.
  [DOI: 10.1038/ncomms13394](http://doi.org/10.1038/ncomms13394)
* M. Dutour Sikirić, O. Delgado-Friedrichs, M. Deza,
  Acta Crystallogr. A 2010, 66, 602.
  [DOI: 10.1107/S0108767310022932](http://doi.org/10.1107/S0108767310022932)
* A. Falenty, T. C. Hansen, W. F. Kuhs, 
  Nature 2014, 516, 231.
  [DOI: 10.1038/nature14014](http://doi.org/10.1038/nature14014)
* Xiaofeng Fan, Dan Bing, Jingyun Zhang, Zexiang Shen, Jer-Lai Kuo,
  Computational Materials Science 49 (2010) S170–S175.
* C. J. Fennell, J. D. Gezelter, 
  J. Chem. Theory Comput. 2005, 1, 662.
  [DOI: 10.1021/ct050005s](http://doi.org/10.1021/ct050005s)
* M. Hirata, T. Yagasaki, M. Matsumoto, H. Tanaka, 
  Langmuir 33, 42, 11561-11569(2017).
  [DOI: 10.1021/acs.langmuir.7b01764](http://doi.org/10.1021/acs.langmuir.7b01764)
* Y. Huang, C. Zhu, L. Wang, X. Cao, Y. Su, X. Jiang, S. Meng, J. Zhao, X. C. Zeng, 
  Science Advances 2016, 2, e1501010.
  [DOI: 10.1126/sciadv.1501010](http://doi.org/10.1126/sciadv.1501010)
* Y. Huang, C. Zhu, L. Wang, J. Zhao, X. C. Zeng,
  Chem. Phys. Lett. 2017, 671, 186.
  [DOI: 10.1016/j.cplett.2017.01.035](http://doi.org/10.1016/j.cplett.2017.01.035)
* G. A. Jeffrey, 
  In Inclusion Compounds; 
  J. L. Atwood, J. E. D. Davies, D. D. MacNicol, Eds.;
  Academic Press: London, 1984, Vol. 1, Chap. 5.
* A. J. Karttunen, T. F. Fässler, M. Linnolahti, T. A. Pakkanen,
  Inorg Chem 2011, 50, 1733.
  [DOI: 10.1021/ic102178d](http://doi.org/10.1021/ic102178d)
* K Koga, GT Gao, H Tanaka, XC Zeng, 
  Nature 412 (6849), 802-805.
  [DOI:10.1038/35090532](http://doi.org/10.1038/35090532)
* V.I. Kosyakov, T.M. Polyanskaya, 
  J. Struct. Chem. 1999, 40, 239.
  [DOI:10.1007/BF02903652](http://doi.org/10.1007/BF02903652)
* WF Kuhs, JL Finney, C Vettier, DV Bliss,
  J. Chem. Phys. 81, 3612–3623 (1998).
  [DOI: 10.1063/1.448109](http://doi.org/10.1063/1.448109)
* Y. Liu, Y. Huang, C. Zhu, H. Li, J. Zhao, L. Wang, L. Ojamäe, J.S. Francisco, X.C. Zeng,
  PNAS 116, 12684-12691 (2019).
  [DOI: 10.1073/pnas.1900739116](https://doi.org/10.1073/pnas.1900739116)
* C. Lobban, J.L. Finney, W. F. Kuhs, 
  Nature 1998, 391, 268.
  [DOI: 10.1038/34622](http://doi.org/10.1038/34622)
* D. Londono, W.F. Kuhs, J.L. Finney,
  Nature 1988, 332, 141.
  [DOI: 10.1038/332141a0](http://doi.org/10.1038/332141a0)
* D. Londono, W.F. Kuhs, J.L. Finney,
  J. Chem. Phys. 98, 4878–4888 (1993).
  [DOI: 10.1063/1.464942](http://doi.org/10.1063/1.464942)
* T. Matsui, M. Hirata, T. Yagasaki, M. Matsumoto, H. Tanaka,
  J. Chem. Phys. 147, 091101 (2017).
  [DOI: /10.1063/1.4994757](http://doi.org/10.1063/1.4994757)
* T. Matsui, T. Yagasaki, M. Matsumoto, H. Tanaka,
  J. Chem. Phys. 150, 041102 (2019).
  [DOI: 10.1063/1.5083021](http://doi.org/10.1063/1.5083021)
* M. Matsumoto, T. Yagasaki, H. Tanaka,
  J. Chem. Phys. 150, 214504 (2019).
  [DOI: 10.1063/1.5096556](https://doi.org/10.1063/1.5096556)
* K. Mochizuki, K. Himoto, M. Matsumoto,
  Phys. Chem. Chem. Phys. 16, 16419–16425 (2014).
  [DOI: 10.1039/c4cp01616e](http://doi.org/10.1039/c4cp01616e)
* N. Mousseau, G. T. Barkema, 
  Curr. Opin. Solid State Mater. Sci. 2001, 5, 497.
  [DOI: 10.1016/S1359-0286(02)00005-0](https://doi.org/10.1016/S1359-0286%2802%2900005-0)
* T. Nakamura, M. Matsumoto, T. Yagasaki, H. Tanaka, 
  J. Phys. Chem. B 2015, 120, 1843.
  [DOI: 10.1021/acs.jpcb.5b09544](http://doi.org/10.1021/acs.jpcb.5b09544)
* J. Russo, F. Romano, H. Tanaka,
  Nat. Mater. 2014, 13, 733.
  [DOI: 10.1038/NMAT3977](http://doi.org/10.1038/NMAT3977)
* C.G. Salzmann, P. Radaelli, A. Hallbrucker, E. Mayer,
  Science 311, 1758–1761 (2006).
  [DOI: 10.1126/science.1123896](http://doi.org/10.1126/science.1123896)
* G. S. Smirnov, V. V. Stegailov,
  J. Phys. Chem. Lett. 2013, 4, 3560.
  [DOI: 10.1021/jz401669d](http://doi.org/10.1021/jz401669d)
* W. L. Vos, L. W. Finger, R. J. Hemley, H. Mao,
  Phys. Rev. Lett. 1993, 71, 3150.
  [DOI: 10.1103/PhysRevLett.71.3150](http://doi.org/10.1103/PhysRevLett.71.3150)
* T. Yagasaki, M. Matsumoto, H. Tanaka,
  J. Phys. Chem. B 122, 7718–7725 (2018).
  [DOI: 10.1021/acs.jpcb.8b04441](http://doi.org/10.1021/acs.jpcb.8b04441)

# Algorithm and how to cite it.

The algorithm to make a depolarized hydrogen-disordered ice is explained in our recent paper:

M. Masakazu, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered
Ice Generator",  J. Comput. Chem. 39, 61-64 (2017). [DOI: 10.1002/jcc.25077](http://doi.org/10.1002/jcc.25077)

    @article{Matsumoto:2017bk,
        author = {Matsumoto, Masakazu and Yagasaki, Takuma and Tanaka, Hideki},
        title = {{GenIce: Hydrogen-Disordered Ice Generator}},
        journal = {Journal of Computational Chemistry},
		volume = {39},
		pages = {61-64},
        year = {2017}
    }

# How to contribute

GenIce is served at the GitHub (https://github.com/vitroid/GenIce/) as an open source software since 2015. Feedbacks, proposals for improvements and extensions, and bug fixes are sincerely appreciated. Developers and test users are also welcome. Please let us know if there are ices that have been published but is not in GenIce.
