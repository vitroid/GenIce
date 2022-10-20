![Logo](https://raw.githubusercontent.com/vitroid/GenIce/develop/logo/genice-v0.png)
# GenIce2

A Swiss army knife to generate hydrogen-disordered ice structures.

version 2.1.5


## New in GenIce2.1

GenIce2-MDAnalysis integration is now available. Try
```shell
% pip install genice2-mdanalysis
% genice2 1h -r 4 4 4 -f "mdanalysis[1h.pdb]"
```
to generate a PDB file.

## Demo

The new GenIce works very well with interactive execution.
[Try instantly](https://colab.research.google.com/github/vitroid/GenIce/blob/main/jupyter.ipynb) on Google Colab.

## Requirements

* networkx>=2.0.dev20160901144005
* cycless
* pairlist>=0.2.12.4
* yaplotlib>=0.1
* numpy
* wheel
* openpyscad
* graphstat
* tilecycles>=0.1.4.3


**Note**: In case you encounter an error complaining "No module named '_ctypes'": Python3.7 and later may require `libffi` for `pairlist` and `tilecycles` modules. Please install `libffi-devel` via the package management system for your system (apt, yum, dnf, brew, etc.)

**Note 2**: There may be compatibility issues when you install GenIce in Apple M1.

1. networkx requires scipy, which cannot be installed simply. See thr workaround at https://stackoverflow.com/questions/65745683/how-to-install-scipy-on-apple-silicon-arm-m1

2. Moreover, scipy requires pythran (I do not know what it is).

```shell
pip3 install pythran
pip3 install cython pybind11
pip3 install --no-binary :all: --no-use-pep517 numpy
brew install openblas gfortran
export OPENBLAS=/opt/homebrew/opt/openblas/lib/
pip3 install --no-binary :all: --no-use-pep517 scipy
```


## Installation
GenIce is registered to [PyPI (Python Package Index)](https://pypi.python.org/pypi/GenIce).
Install with pip3.

    pip3 install genice2

## Uninstallation

    pip3 uninstall genice2

## Usage

    usage: genice2 [-h] [--version] [--rep REP REP REP]
                   [--shift SHIFT SHIFT SHIFT] [--dens DENS] [--add_noise percent]
                   [--seed SEED] [--format name] [--water model] [--guest 14=me]
                   [--Guest 13=me] [--Group 13=bu-:0] [--anion 3=Cl]
                   [--cation 3=Na] [--depol DEPOL] [--asis] [--debug] [--quiet]
                   [--assess_cages]
                   Type

    GenIce is a swiss army knife to generate hydrogen-disordered ice structures.
    (version 2.1.5)

    positional arguments:
      Type                  Crystal type (1c, 1h, etc. See
                            https://github.com/vitroid/GenIce for available ice
                            structures.)
                            If you want to analyze your own structures, please try
                            analice tool.


                            [Available lattice structures]

                            1. Lattice structures served with GenIce

                            0, ice0         Metastable ice "0".
                            11, XI, ice11   A candidate for an antiferroelectric
                                            Ice XI #19.
                            115_2_114, 12_1_11, 144_2_7301, 151_2_4949650,
                                            153_2_155471, 176_2_5256, 207_1_4435,
                                            2_2_623457, ACO, CS4, DDR, IWV, LTA,
                                            MAR, NON, PCOD8007225, PCOD8036144,
                                            PCOD8204698, PCOD8301974, PCOD8321499,
                                            PCOD8324623, SGT, SOD, engel01,
                                            engel03, engel04, engel17, engel20,
                                            engel23, engel24, engel26, engel29,
                                            engel30, engel31, engel34, sVII
                                            Hypothetical zeolitic ice
                            11i             Sixteen candidates for Ice XI.
                            12, XII, ice12  Metastable high-pressure ice XII.
                            13, XIII, ice13 Ice XIII, a hydrogen-ordered
                                            counterpart of ice V.
                            16, CS2, MTN, XVI, ice16, sII   Ultralow-density Ice
                                            XVI.
                            17, XVII, ice17 Ultralow-density Ice XVII.
                            1c, Ic, ice1c   Cubic type of ice I.
                            1h, Ih, ice1h   Most popular Ice I (hexagonal)
                            2, II, ice2     Hydrogen-ordered ice II.
                            2D3             Trilayer honeycomb ice.
                            2d, ice2d, ice2rect             A hydrogen-disordered
                                            counterpart of ice II.
                            3, III, ice3    Ice III.
                            4, IV, ice4     Ice IV.
                            4R              Orthogonalized ice IV.
                            5, V, ice5      Monoclinic ice V (testing).
                            5R              Ice V with orthogonal unit cell.
                                            (testing)
                            6, VI, ice6     Conventional high-pressure ice VI.
                            6h              Half lattice of ice VI.
                            7, VII, ice7    Conventional high-pressure ice VII.
                            8, VIII, ice8   Ice VIII, a hydrogen-ordered
                                            counterpart of ice VII.
                            9, IX, ice9     Ice IX, a hydrogen-ordered counterpart
                                            of ice III.
                            A, iceA         Hypothetical ice A.
                            A15, Struct33   Cubic Structure I of clathrate hydrate.
                            B, iceB         Hypothetical ice B.
                            BSV, engel05    Hypothetical zeolitic ice of the gyroid
                                             structure.
                            C14, C15, C36, FK6layers, FK9layers, HS2, Hcomp,
                                            Struct01, Struct03, Struct04, Struct05,
                                             Struct06, Struct07, Struct08,
                                            Struct09, Struct10, Struct11, Struct12,
                                             Struct13, Struct14, Struct15,
                                            Struct16, Struct17, Struct18, Struct19,
                                             Struct20, Struct21, Struct22,
                                            Struct23, Struct24, Struct25, Struct26,
                                             Struct27, Struct28, Struct29,
                                            Struct30, Struct31, Struct32, Struct34,
                                             Struct35, Struct36, Struct37,
                                            Struct38, Struct39, Struct40, Struct41,
                                             Struct42, Struct43, Struct44,
                                            Struct45, Struct46, Struct47, Struct48,
                                             Struct49, Struct50, Struct51,
                                            Struct52, Struct53, Struct54, Struct55,
                                             Struct56, Struct57, Struct58,
                                            Struct59, Struct60, Struct61, Struct62,
                                             Struct63, Struct64, Struct65,
                                            Struct66, Struct67, Struct68, Struct69,
                                             Struct70, Struct71, Struct72,
                                            Struct73, Struct74, Struct75, Struct76,
                                             Struct77, Struct78, Struct79,
                                            Struct80, Struct81, Struct82, Struct83,
                                             Struct84, Z, delta, mu, psigma, sV,
                                            sigma, zra-d     A space fullerene.
                            CRN1, CRN2, CRN3                A continuous random
                                            network of Sillium.
                            CS1, MEP, sI    Clathrate hydrates sI.
                            DOH, HS3, sH    Clathrate type H.
                            EMT             Hypothetical ice with a large cavity.
                            FAU             Hypothetical ice at negative pressure
                                            ice 'sIV'.
                            RHO             Hypothetical ice at negative pressure
                                            ice 'sIII'.
                            Struct02        A space fullerene. (I phase?)
                            T               Hypothetical clathrate type T.
                            bilayer         A Bilayer Honeycomb Ice Phase in
                                            Hydrophobic Nanopores.
                            c0te            Filled ice C0 by Teeratchanan
                                            (Hydrogen-disordered.) (Positions of
                                            guests are supplied.)
                            c1te            Hydrogen-ordered hydrogen hydrate C1 by
                                             Teeratchanan. (Positions of guests are
                                             supplied.)
                            c2te            Filled ice C2 (cubic ice) by
                                            Teeratchanan (Hydrogen disordered).
                                            (Positions of guests are supplied.)
                            eleven          Ice XI w/ stacking faults.
                            i               Hypothetical ice "i".
                            ice1hte         Filled ice Ih by Teeratchanan (Hydrogen
                                             disordered). (Positions of guests are
                                            supplied.)
                            iceR            Hypothetical ice R.
                            iceT            Hypothetical ice T.
                            iceT2           Hypothetical ice T2.
                            one             Ice I w/ stacking faults.
                            oprism          Hydrogen-ordered ice nanotubes.
                            sTprime         Filled ice sT'.
                            xFAU            Aeroice xFAU.
                            xdtc            A porous ice with cylindrical channels.
                            ----
                            (Undocumented) 1h_unit HS1 TS1 dtc ice1h_unit sIII sIV


                            2. Lattice structures served by external plugins

                            (None)
                            ----
                            (Undocumented) cif  zeolite


                            3. Lattice structures served locally

                            (None)
                            ----



    optional arguments:
      -h, --help            show this help message and exit
      --version, -V         show program's version number and exit
      --rep REP REP REP, -r REP REP REP
                            Repeat the unit cell along a, b, and c axes. [1,1,1]
      --shift SHIFT SHIFT SHIFT, -S SHIFT SHIFT SHIFT
                            Shift the unit cell along a, b, and c axes. (0.5==half
                            cell) [0,0,0]
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
                            raw             Raw data. (For use with Jupyter)
                            rings           Show rings in Yaplot.
                            y, yaplot       Yaplot.
                            ----
                            (Undocumented) c cif cif2 com null r rcom towhee xyz


                            2. Formatters served by external plugins

                            (None)
                            ----
                            (Undocumented) _RDF  png  svg


                            3. Formatters served locally

                            (None)
                            ----


      --water model, -w model
                            Specify the water model. [tip3p]


                            [Available molecules]

                            1. Molecules served with GenIce

                            3site, tip3p    A typical 3-site model.
                            4site, tip4p    A typical 4-site model.
                            5site, tip5p    A typical 5-site model.
                            6site, NvdE     A 6-site water model.
                            7site           A seven-site water model.
                            physical_water  Physical model of water; Oxygen atom is
                                             on the lattice point.
                            ----
                            (Undocumented) ice spce


                            2. Molecules served by external plugins

                            (None)
                            ----


                            3. Molecules served locally

                            (None)
                            ----


      --guest 14=me, -g 14=me
                            Specify guest(s) in the cage type. (D=empty,
                            T=co2*0.5+me*0.3, etc.)


                            [Available molecules]

                            1. Molecules served with GenIce

                            H2              Hydrogen molecule.
                            ch4             An all-atom methane model.
                            et              A united-atom ethane model.
                            me              A united-atom methane model.
                            mol             Loader for MOL files (generated by
                                            MolView.org), e.g. mol[THF.mol].
                            thf             An all-atom tetrahydrofuran (THF)
                                            model.
                            uathf           A united-atom five-site tetrahydrofuran
                                             (THF) model.
                            ----
                            (Undocumented) co2 empty g12 g14 g15 g16 one uathf6


                            2. Molecules served by external plugins

                            (None)
                            ----


                            3. Molecules served locally

                            (None)
                            ----


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
      --depol DEPOL         Depolarization. (strict, optimal, or none) ["strict"]
      --asis                Assumes all given HB pairs to be fixed. No shuffle and
                            no depolarization.
      --debug, -D           Output debugging info.
      --quiet, -q           Do not output progress messages.
      --assess_cages, -A    Assess the locations of cages based on the HB network
                            topology. Note: it may fail when the unit cell is too
                            small.


Use `./genice.x` instead of `genice2` if you want to use it inside the source tree.

## Examples

* To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
.gro format:

        genice2 --water tip4p --rep 3 3 3  4 > ice4.gro

* To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing THF (united atom with a dummy site) in the large cage in GROMACS
.gro format:

        genice2 -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro


## Basics

The program generates various ice lattice with proton disorder and without defect.  The total dipole moment is always set to zero (except in the case you specify `--depol` option).  The minimal structure (with --rep 1 1 1 option) is not always the unit cell of the lattice because it is difficult to deal with the hydrogen bond network topology of tiny lattice under periodic boundary condition.  Note that the generated structure is not optimal according to the potential energy.

* To get a large repetition of ice Ih in XYZ format,

        genice2 --rep 8 8 8 1h --format xyz > 1hx888.xyz

* To get a ice V lattice of different hydrogen order in CIF format, use `-s` option to specify the random seed.

        genice2 5 -s 1024 --format cif > 5-1024.cif

* To obtain an ice VI lattice with different density and with TIP4P water model in gromacs format, use `--dens x` option to specify the density in g cm<sup>-3</sup>.

        genice2 6 --dens 1.00 --format g --water tip4p > 6d1.00.gro

GenIce is a modular program; it reads a unit cell data from a lattice plugin defined in the lattices folder, put water and guest molecules using a molecule plugin defined in the molecules/ folder, and output in various formats using a format plugin defined in the formats/ folder. You can write your plugins to extend GenIce. Some plugins also accept options.

## Clathrate hydrates

For clathrate hydrates, you can prepare the lattice with cages partially occupied by various guest molecules.

* To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS .gro format: (60% of small cages are filled with co2 and 40% are methane)

        genice2 -g 12=co2*0.6+me*0.4 -g 14=co2 --water tip4p CS1 > cs1.gro

* To make a CS2 clathrate hydrate structure of TIP5P water containing THF molecules in the large cage, while only one cage is filled with methane molecule, first, just run `genice2` without guest specifications:

        genice2 CS2 > CS2.gro

    The list of cages will be output as follows:

        INFO   Cage types: ['12', '16']
        INFO   Cage type 12: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183}
        INFO   Cage type 16: {136, 137, 138, 139, 140, 141, 142, 143, 16, 17, 18, 19, 20, 21, 22, 23, 160, 161, 162, 163, 164, 165, 166, 167, 40, 41, 42, 43, 44, 45, 46, 47, 184, 185, 186, 187, 188, 189, 190, 191, 64, 65, 66, 67, 68, 69, 70, 71, 88, 89, 90, 91, 92, 93, 94, 95, 112, 113, 114, 115, 116, 117, 118, 119}

    This indicates that there are two types of cages named `12` and `16`.  Fill the `16` cages with THF and put a methane molecule in the `0`th cage of type `12` as follows:

        genice2 CS2 -g 16=uathf -G 0=me > CS2.gro

Although only a few kinds of guest molecules are preset, you can easily prepare new guest molecules as a module. Here is an example of the ethylene oxide molecule.

```python
import numpy as np
import genice2.molecules

class Molecule(genice2.molecules.Molecule):
    def __init__(self):
        # United-atom EO model with a dummy site
        LOC = 0.1436 # nm
        LCC = 0.1472 # nm

        Y = (LOC**2 - (LCC/2)**2)**0.5

        self.sites_ = np.array([[ 0.,    0., 0. ],
                               [-LCC/2, Y,  0. ],
                               [+LCC/2, Y,  0. ],])

        mass = np.array([16,14,14])
        # center of mass
        CoM = mass @ self.sites / np.sum(mass)
        self.sites_ -= CoM

        self.atoms_  = ["O","C","C"]
        self.labels_ = ["Oe","Ce","Ce"]
        self.name_   = "EO"
```

Write the code in eo.py. Make a folder named `molecules` in the current working directory and put it in.

*Note*: multiple occupancy is not implemented. If it is required, make a module of a virtual molecule that contains multiple molecules.

## Doping ions

Small ions may replace the host molecules.  In that case, you can use `-a` and `-c` options to replace the specified water molecules with anions and cations.

The following example replaces the `0`th water molecule (in the replicated lattice) with Na cation and `1`st water molecule with Cl anion.  The hydrogen bonds around the ions are organized appropriately.

    genice2 CS2 --depol=optimal -c 0=Na -a 1=Cl > CS2.gro

*Note 1*: The numbers of cations and anions must be the same.  Otherwise, the ice rule is never satisfied and the program does not stop.

*Note 2*: The option `--depol=optimal` is also required because it is impossible to completely depolarize the structure containing ions.

*Note 3*: Protonic defects (H<sub>3</sub>O<sup>+</sup> and OH<sup>-</sup>) are not yet implemented.

## Semiclathrate hydrates

### Placement of a tetrabutylammonium ion (testing)

Let us assume that the id of the water molecule to be replaced by nitrogen of the TBA as zero.  Place the nitrogen as a cation and also replace the water 2 by the counter-ion Br.

    genice2 HS1 -c 0=N -a 2=Br --depol=optimal > HS1.gro

Then you will see the following info.

```
INFO   Hints:
INFO     Cage types: ['12', '14', '15']
INFO     Cage type 12: {0, 1, 2, 3, 4, 5}
INFO     Cage type 14: {8, 9, 6, 7}
INFO     Cage type 15: {10, 11, 12, 13}
...
INFO Stage7: Arrange guest atoms.
INFO     Cages adjacent to dopant 2: {0, 9, 2, 13}
INFO     Cages adjacent to dopant 0: {0, 9, 2, 7}
```

It indicates that the nitrogen is surrounded by cages with ids 0, 9, 2, and 7.  Types for these cages can also be found in the info.  Then, we put the Bu- group (minus does not mean ions) in these cages adjacent to dopant 0.

```shell
genice2 HS1 -c 0=N -a 2=Br -H 0=Bu-:0 -H 9=Bu-:0 -H 2=Bu-:0 -H 7=Bu-:0 --depol=optimal > HS1.gro
```

Here the option `-H` specifies the group by `-H (cage id)=(group name):(root)`, and the root is the nitrogen that is specified by `-c` (cation) option.


### Placement of TBAB in the lattice module

*Under preparation*

It is more convenient if the lattice of the semiclathrate hydrate contains molecular ions in the appropriate locations in advance.  Here we explain the way to make the special module for semclathrates.

## Output formats


Name |Application | extension | water | solute | HB | remarks
-------|------------|-----------|----------|---------|-----|---
`cif`, `cif2` |CIF         | `.cif`      | Atomic positions | Atomic positions | none |Experimental
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

By installing the [`genice2-mdanalysis`](https://github.com/vitroid/genice-mdanalysis) package separately, you can generate files in many formats for a  large number of molecular dynamics package softwares. E.g.

```shell
% pip install genice2-mdanalysis
% genice2 1c -f 'mdanalysis[1c.pdb]'
% genice2 1h -f 'mdanalysis[1h.xtc]'
```

All the supported file types are listed in the [MDAnalysis web page](https://docs.mdanalysis.org/stable/documentation_pages/coordinates/init.html#supported-coordinate-formats).

You can prepare your file formats. Create a folder named `formats` in the current working directory and put the plugins in it.

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

Symbol | <div style="width:300px">Description</div>
-------|------------
0, ice0 | Metastable ice "0". [Russo 2014]
11, XI, ice11 | A candidate for an antiferroelectric Ice XI #19. [Jackson 1997, Fan 2010]
115_2_114, 12_1_11, 144_2_7301, 151_2_4949650, 153_2_155471, 176_2_5256, 207_1_4435, 2_2_623457, ACO, CS4, DDR, IWV, LTA, MAR, NON, PCOD8007225, PCOD8036144, PCOD8204698, PCOD8301974, PCOD8321499, PCOD8324623, SGT, SOD, engel01, engel03, engel04, engel17, engel20, engel23, engel24, engel26, engel29, engel30, engel31, engel34, sVII | Hypothetical zeolitic ice [Jeffrey 1984, Kosyakov 1999, Engel 2018, IZA Database]
11i | Sixteen candidates for Ice XI. [Hirsch 2004]
12, XII, ice12 | Metastable high-pressure ice XII. [Lobban 1998, Koza 2000]
13, XIII, ice13 | Ice XIII, a hydrogen-ordered counterpart of ice V. [Salzmann 2006]
16, CS2, MTN, XVI, ice16, sII | Ultralow-density Ice XVI. [Jeffrey 1984, Kosyakov 1999, Sikiric 2010, Falenty 2014, IZA Database]
17, XVII, ice17 | Ultralow-density Ice XVII. [Smirnov 2013, Strobel 2016, Rosso 2016]
1c, Ic, ice1c | Cubic type of ice I. [Vos 1993]
1h, Ih, ice1h | Most popular Ice I (hexagonal)
2, II, ice2 | Hydrogen-ordered ice II. [Kamb 1964, Londono 1988, Kamb 2003]
2D3 | Trilayer honeycomb ice.
2d, ice2d, ice2rect | A hydrogen-disordered counterpart of ice II. [Nakamura 2015]
3, III, ice3 | Ice III. [Petrenko 1999]
4, IV, ice4 | Ice IV. [Avogadro]
4R | Orthogonalized ice IV. [Avogadro]
5, V, ice5 | Monoclinic ice V (testing).
5R | Ice V with orthogonal unit cell. (testing)
6, VI, ice6 | Conventional high-pressure ice VI. [Petrenko 1999]
6h | Half lattice of ice VI.
7, VII, ice7 | Conventional high-pressure ice VII.
8, VIII, ice8 | Ice VIII, a hydrogen-ordered counterpart of ice VII. [Kuhs 1998]
9, IX, ice9 | Ice IX, a hydrogen-ordered counterpart of ice III. [Londono 1993]
A, iceA | Hypothetical ice A. [Baez 1998]
A15, Struct33 | Cubic Structure I of clathrate hydrate. [Sikiric 2010]
B, iceB | Hypothetical ice B. [Baez 1998]
BSV, engel05 | Hypothetical zeolitic ice of the gyroid structure. [Engel 2018, IZA Database]
C14, C15, C36, FK6layers, FK9layers, HS2, Hcomp, Struct01, Struct03, Struct04, Struct05, Struct06, Struct07, Struct08, Struct09, Struct10, Struct11, Struct12, Struct13, Struct14, Struct15, Struct16, Struct17, Struct18, Struct19, Struct20, Struct21, Struct22, Struct23, Struct24, Struct25, Struct26, Struct27, Struct28, Struct29, Struct30, Struct31, Struct32, Struct34, Struct35, Struct36, Struct37, Struct38, Struct39, Struct40, Struct41, Struct42, Struct43, Struct44, Struct45, Struct46, Struct47, Struct48, Struct49, Struct50, Struct51, Struct52, Struct53, Struct54, Struct55, Struct56, Struct57, Struct58, Struct59, Struct60, Struct61, Struct62, Struct63, Struct64, Struct65, Struct66, Struct67, Struct68, Struct69, Struct70, Struct71, Struct72, Struct73, Struct74, Struct75, Struct76, Struct77, Struct78, Struct79, Struct80, Struct81, Struct82, Struct83, Struct84, Z, delta, mu, psigma, sV, sigma, zra-d | A space fullerene. [Sikiric 2010]
CRN1, CRN2, CRN3 | A continuous random network of Sillium. [Mousseau 2001]
CS1, MEP, sI | Clathrate hydrates sI. [Frank 1959, Jeffrey 1984, Kosyakov 1999, IZA Database]
DOH, HS3, sH | Clathrate type H.
EMT | Hypothetical ice with a large cavity. [Liu 2019, IZA Database]
FAU | Hypothetical ice at negative pressure ice 'sIV'. [Huang 2017, IZA Database]
RHO | Hypothetical ice at negative pressure ice 'sIII'. [Huang 2016, IZA Database]
Struct02 | A space fullerene. (I phase?) [Sikiric 2010]
T | Hypothetical clathrate type T. [Sikiric 2010, Karttunen 2011]
bilayer | A Bilayer Honeycomb Ice Phase in Hydrophobic Nanopores. [Koga 1997]
c0te | Filled ice C0 by Teeratchanan (Hydrogen-disordered.) (Positions of guests are supplied.) [Teeratchanan 2015]
c1te | Hydrogen-ordered hydrogen hydrate C1 by Teeratchanan. (Positions of guests are supplied.) [Teeratchanan 2015]
c2te | Filled ice C2 (cubic ice) by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.) [Teeratchanan 2015]
eleven | Ice XI w/ stacking faults.
i | Hypothetical ice "i". [Fennell 2005]
ice1hte | Filled ice Ih by Teeratchanan (Hydrogen disordered). (Positions of guests are supplied.) [Teeratchanan 2015]
iceR | Hypothetical ice R. [Maynard-Casely 2010, Mochizuki 2014]
iceT | Hypothetical ice T. [Hirata 2017]
iceT2 | Hypothetical ice T2. [Yagasaki 2018]
one | Ice I w/ stacking faults.
oprism | Hydrogen-ordered ice nanotubes. [Koga 2001]
sTprime | Filled ice sT'. [Smirnov 2013]
xFAU | Aeroice xFAU. [Matsui 2017]
xdtc | A porous ice with cylindrical channels. [Matsumoto 2021]
1h_unit, HS1, TS1, dtc, ice1h_unit, sIII, sIV | (Undocumented)


Ice names with double quotations are not experimentally verified.

You can prepare your own ice structures. Create a folder named `lattices` in the current working directory and put the plugins in it.

[cif2ice](https://github.com/vitroid/cif2ice) is a tool to retrieve a
cif file of zeolite from [IZA structure database](http://www.iza-structure.org/databases) and prepare a lattice
module in the path above.

Note: Different names are given in different nomenclature.

CH/FI|CH  |ice|FK |Zeo|Foam|
-----|----|---|---|---|----|
sI   |CS1 |-  |A15|MEP|Weaire-Phelan|
sII  |CS2 |16 |C15|MTN|    |
sIII |TS1 |-  |sigma|-  |    |
sIV  |HS1 |-  |Z  |-  |    |
sV   |HS2 |-  |*  |-  |    |
sVII |CS4 |-  |*  |SOD|Kelvin|
sH   |HS3 |-  |*  |DOH|    |
C0   |-   |17 |*  |-  |    |
C1   |-   |2  |*  |-  |    |
C2   |-   |1c |*  |-  |    |

FI: Filled ices; CH: Clathrate hydrates; FK:Frank-Kasper duals; Zeo: Zeolites; Foam: foam crystals (Weaire 1994).

-: No correspondence; *: Non-FK types.

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.
## Water models
A water model can be chosen with `--water` option.

symbol   | type
---------|-----
3site, tip3p | A typical 3-site model.
4site, tip4p | A typical 4-site model. [Jorgensen 1983, Jorgensen 1985]
5site, tip5p | A typical 5-site model.
6site, NvdE | A 6-site water model. [Nada 2003]
7site | A seven-site water model. [Zhao 2019]
physical_water | Physical model of water; Oxygen atom is on the lattice point. [Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem Phys, 79, 926 (1983).]
ice, spce | (Undocumented)


## Guest molecules

symbol | type
-------|---------
H2 | Hydrogen molecule. [https://www.britannica.com/science/hydrogen]
ch4 | An all-atom methane model.
et | A united-atom ethane model.
me | A united-atom methane model.
mol | Loader for MOL files (generated by MolView.org), e.g. mol[THF.mol].
thf | An all-atom tetrahydrofuran (THF) model.
uathf | A united-atom five-site tetrahydrofuran (THF) model.
co2, empty, g12, g14, g15, g16, one, uathf6 | (Undocumented)


You can prepare your own guest molecules.  Create a folder named `molecules` in the current working directory and put the plugins in it.

# Extra plugins

Some extra plugins are available via python package index using pip command.

For example, you can install RDF plugin by the following command,

    % pip install genice2-rdf

And use it as an output format to get the radial distribution functions.

    % genice2 TS1 -f _RDF > TS1.rdf.txt



## Output and analysis plugins

Analysis plugin is a kind of output plugin (specified with -f option).

| pip name | GenIce2 option | Description | output format | requirements |
|----------|-------|-------------|---------------|--------------|
|[`genice2-cage`](https://github.com/vitroid/genice-cage)|`-f _cage`| Detect cages and quasi-polyhedra (vitrites). | text, json, gromacs | `cycless` |
|[`genice2-rdf`](https://github.com/vitroid/genice-rdf)|`-f _RDF`| Radial distribution functions. | text |  |
|[`genice2-svg`](https://github.com/vitroid/genice-svg)|`-f svg`<br />`-f png` | 2D graphics in SVG format.<br /> ... in PNG format.| SVG<br />PNG | `svgwrite` |
|[`genice2-twist`](https://github.com/vitroid/genice-twist)|`-f twist`| Calculate the twist order parameter (and visualize) [Matsumoto 2019]| text (@BTWC)<br />SVG<br />PNG <br />yaplot | `twist-op`, `genice2-svg` |
|[`genice2-mdanalysis`](https://github.com/vitroid/genice-mdanalysis)|`-f mdanalysis`| Output the atoms in various file formats that are served by [MDAnalysis](https://github.com/MDAnalysis/mdanalysis).| text, binary | `mdanalysis` |

## Input plugins

Input plugins (a.k.a. lattice plugins) construct a crystal structure on demand.

| pip name   | GenIce2 usage    | Description  |requirements |
|------------|-----------------|--------------|-------------|
|[`genice2-cif`](https://github.com/vitroid/genice-cif)| `genice2 cif[ITT.cif]`<br /> `genice2 zeolite[ITT]`| Read a local CIF file as an ice structure.<br />Read a structure from Zeolite DB. | `cif2ice` |

# Changes from GenIce1

## Novel algorithm to make a structure obeying the ice rules in Stage 3.

- We have devised a completely new algorithm for orienting water molecules so that they obey ice rules. This algorithm can be applied only to defect-free ice. The algorithm runs in the following steps.
  1. First, based on the distances between neighboring molecules, the structure of the hydrogen-bond network is represented by a 4-connected undirected graph.
  2. The undirected graph is then randomly tiled with cycles. That is, we draw many cycles in the network so that all edges belong to only one of the cycles. It is always possible to a 4-connected regular graph.
  3. By directing each cycle, we can immediately obtain a directed graph that satisfies the ice rule. We can choose two orientations for each cycle so that the total polarization of the entire system is as small as possible.
  4. In rare cases, complete depolarization may not be possible. In such cases, it is depolarized in Stage 4.

## Faster, faster, faster.

Combinations of the new algorithm and other improvements in coding, the processing time of GenIce2 is about five times faster than that of GenIce1.

## Core algorithm is separated.

The core part of the new algorithm is separated as the TileCycles package.

## Colaboratory-ready!

Now GenIce2 works on the [Google Colaboratory!](https://colab.research.google.com/github/vitroid/GenIce/blob/genice2/jupyter.ipynb)

## New ices

Many new ice structures are added.

## Integration with MDAnalysis

GenIce2 is now integrated with MDAnalysis.

# References

* [Avogadro] Avogadro https://github.com/cryos/avogadro/blob/master/crystals/ice/H2O-Ice-IV.cif
* [Baez 1998]
BÁEZ, Luis A. and CLANCY, Paulette, 1995, Phase equilibria in extended simple point charge ice‐water systems. The Journal of Chemical Physics [online]. 8 December 1995. Vol. 103, no. 22, p. 9744–9755. DOI 10.1063/1.469938. Available from: http://dx.doi.org/10.1063/1.469938

* [Cockayne 1995]
COCKAYNE, Eric, 1995. Dense quasiperiodic decagonal disc packing. Physical Review B [online]. 1 June 1995. Vol. 51, no. 21, p. 14958–14961. DOI 10.1103/physrevb.51.14958. Available from: http://dx.doi.org/10.1103/PhysRevB.51.14958

* [Engel 2018]
ENGEL, Edgar A., ANELLI, Andrea, CERIOTTI, Michele, PICKARD, Chris J. and NEEDS, Richard J., 2018, Mapping uncharted territory in ice from zeolite networks to ice structures. Nature Communications [online]. 5 June 2018. Vol. 9, no. 1. DOI 10.1038/s41467-018-04618-6. Available from: http://dx.doi.org/10.1038/s41467-018-04618-6

* [Frank 1959]
FRANK, F. C. and KASPER, J. S., 1959, Complex alloy structures regarded as sphere packings. II. Analysis and classification of representative structures. Acta Crystallographica [online]. 10 July 1959. Vol. 12, no. 7, p. 483–499. DOI 10.1107/s0365110x59001499. Available from: http://dx.doi.org/10.1107/S0365110X59001499

* [Falenty 2014]
FALENTY, Andrzej, HANSEN, Thomas C. and KUHS, Werner F., 2014, Formation and properties of ice XVI obtained by emptying a type sII clathrate hydrate. Nature [online]. December 2014. Vol. 516, no. 7530, p. 231–233. DOI 10.1038/nature14014. Available from: http://dx.doi.org/10.1038/nature14014

* [Fan 2010] Xiaofeng Fan, Dan Bing, Jingyun Zhang, Zexiang Shen, Jer-Lai Kuo,  Computational Materials Science 49 (2010) S170–S175.
* [Fennell 2005]
FENNELL, Christopher J. and GEZELTER, J. Daniel, 2005, Computational Free Energy Studies of a New Ice Polymorph Which Exhibits Greater Stability than Ice Ih. Journal of Chemical Theory and Computation [online]. 30 April 2005. Vol. 1, no. 4, p. 662–667. DOI 10.1021/ct050005s. Available from: http://dx.doi.org/10.1021/ct050005s

* [Hirata 2017]
HIRATA, Masanori, YAGASAKI, Takuma, MATSUMOTO, Masakazu and TANAKA, Hideki, 2017, Phase Diagram of TIP4P/2005 Water at High Pressure. Langmuir [online]. 24 August 2017. Vol. 33, no. 42, p. 11561–11569. DOI 10.1021/acs.langmuir.7b01764. Available from: http://dx.doi.org/10.1021/acs.langmuir.7b01764

* [Hirsch 2004] Hirsch, T. K. & Ojamäe, L. Quantum-Chemical and Force-Field Investigations of Ice Ih: Computation of Proton-Ordered Structures and Prediction of Their Lattice Energies. J. Phys. Chem. B 108, 15856-15864 (2004)
* [Huang 2016]
HUANG, Yingying, ZHU, Chongqin, WANG, Lu, CAO, Xiaoxiao, SU, Yan, JIANG, Xue, MENG, Sheng, ZHAO, Jijun and ZENG, Xiao Cheng, 2016, A new phase diagram of water under negative pressure: The rise of the lowest-density clathrate s-III. Science Advances [online]. February 2016. Vol. 2, no. 2, p. e1501010. DOI 10.1126/sciadv.1501010. Available from: http://dx.doi.org/10.1126/sciadv.1501010

* [Huang 2017]
HUANG, Yingying, ZHU, Chongqin, WANG, Lu, ZHAO, Jijun and ZENG, Xiao Cheng, 2017, Prediction of a new ice clathrate with record low density: A potential candidate as ice XIX in guest-free form. Chemical Physics Letters [online]. March 2017. Vol. 671, p. 186–191. DOI 10.1016/j.cplett.2017.01.035. Available from: http://dx.doi.org/10.1016/j.cplett.2017.01.035

* [IZA Database] http://www.iza-structure.org/databases/
* [Jackson 1997] Jackson, S. M., V. M. Nield, R. W. Whitworth, M. Oguro, and C. C. Wilson, 1997, ‘‘Single-crystal neutron diffraction studies of the structure of ice XI,’’ J. Phys. Chem. B 101, 6142.
* [Jeffrey 1984] G. A. Jeffrey, In Inclusion Compounds; J. L. Atwood, J. E. D. Davies, D. D. MacNicol, Eds.; Academic Press: London, 1984, Vol. 1, Chap. 5.
* [Jorgensen 1983] W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, and M. L. Klein, Comparison of simple potential functions for simulating liquid water, J. Chem. Phys. 79 (1983) 926-935.
* [Jorgensen 1985] W. L. Jorgensen and J. D. Madura, Temperature and size dependence for monte carlo simulations of TIP4P water, Mol. Phys. 56 (1985) 1381-1392.
* [Kamb 1964] Kamb, B.IUCr. Ice. II. A proton-ordered form of ice. Acta Cryst 17, 1437–1449 (1964).
* [Kamb 2003] Kamb, B., Hamilton, W. C., LaPlaca, S. J. & Prakash, A. Ordered Proton Configuration in Ice II, from Single‐Crystal Neutron Diffraction. J. Chem. Phys. 55, 1934–1945 (2003).
* [Karttunen 2011]
KARTTUNEN, Antti J., FÄSSLER, Thomas F., LINNOLAHTI, Mikko and PAKKANEN, Tapani A., 2011, Structural Principles of Semiconducting Group 14 Clathrate Frameworks. Inorganic Chemistry [online]. 7 March 2011. Vol. 50, no. 5, p. 1733–1742. DOI 10.1021/ic102178d. Available from: http://dx.doi.org/10.1021/ic102178d

* [Koga 1997]
KOGA, Kenichiro, ZENG, X. C. and TANAKA, Hideki, 1997. Freezing of Confined Water: A Bilayer Ice Phase in Hydrophobic Nanopores. Physical Review Letters [online]. 29 December 1997. Vol. 79, no. 26, p. 5262–5265. DOI 10.1103/physrevlett.79.5262. Available from: http://dx.doi.org/10.1103/PhysRevLett.79.5262

* [Koga 2001]
KOGA, Kenichiro, GAO, G. T., TANAKA, Hideki and ZENG, X. C., 2001, Formation of ordered ice nanotubes inside carbon nanotubes. Nature [online]. August 2001. Vol. 412, no. 6849, p. 802–805. DOI 10.1038/35090532. Available from: http://dx.doi.org/10.1038/35090532

* [Kosyakov 1999]
KOSYAKOV, V. I. and POLYANSKAYA, T. M., 1999, Using structural data for estimating the stability of water networks in clathrate and semiclathrate hydrates. Journal of Structural Chemistry [online]. March 1999. Vol. 40, no. 2, p. 239–245. DOI 10.1007/bf02903652. Available from: http://dx.doi.org/10.1007/BF02903652

* [Koza 2000]
KOZA, M. M., SCHOBER, H., HANSEN, T., TÖLLE, A. and FUJARA, F., 2000, Ice XII in Its Second Regime of Metastability. Physical Review Letters [online]. 1 May 2000. Vol. 84, no. 18, p. 4112–4115. DOI 10.1103/physrevlett.84.4112. Available from: http://dx.doi.org/10.1103/PhysRevLett.84.4112

* [Kuhs 1998]
KUHS, W. F., FINNEY, J. L., VETTIER, C. and BLISS, D. V., 1984, Structure and hydrogen ordering in ices VI, VII, and VIII by neutron powder diffraction. The Journal of Chemical Physics [online]. 15 October 1984. Vol. 81, no. 8, p. 3612–3623. DOI 10.1063/1.448109. Available from: http://dx.doi.org/10.1063/1.448109

* [Liu 2019]
LIU, Yuan, HUANG, Yingying, ZHU, Chongqin, LI, Hui, ZHAO, Jijun, WANG, Lu, OJAMÄE, Lars, FRANCISCO, Joseph S. and ZENG, Xiao Cheng, 2019, An ultralow-density porous ice with the largest internal cavity identified in the water phase diagram. Proceedings of the National Academy of Sciences [online]. 10 June 2019. Vol. 116, no. 26, p. 12684–12691. DOI 10.1073/pnas.1900739116. Available from: http://dx.doi.org/10.1073/pnas.1900739116

* [Lobban 1998]
LOBBAN, C., FINNEY, J. L. and KUHS, W. F., 1998, The structure of a new phase of ice. Nature [online]. January 1998. Vol. 391, no. 6664, p. 268–270. DOI 10.1038/34622. Available from: http://dx.doi.org/10.1038/34622

* [Londono 1988]
LONDONO, D., KUHS, W. F. and FINNEY, J. L., 1988, Enclathration of helium in ice II: the first helium hydrate. Nature [online]. March 1988. Vol. 332, no. 6160, p. 141–142. DOI 10.1038/332141a0. Available from: http://dx.doi.org/10.1038/332141a0

* [Londono 1993]
LONDONO, J. D., KUHS, W. F. and FINNEY, J. L., 1993, Neutron diffraction studies of ices III and IX on under‐pressure and recovered samples. The Journal of Chemical Physics [online]. 15 March 1993. Vol. 98, no. 6, p. 4878–4888. DOI 10.1063/1.464942. Available from: http://dx.doi.org/10.1063/1.464942

* [Matsui 2017]
MATSUI, Takahiro, HIRATA, Masanori, YAGASAKI, Takuma, MATSUMOTO, Masakazu and TANAKA, Hideki, 2017, Communication: Hypothetical ultralow-density ice polymorphs. The Journal of Chemical Physics [online]. 7 September 2017. Vol. 147, no. 9, p. 091101. DOI 10.1063/1.4994757. Available from: http://dx.doi.org/10.1063/1.4994757

* [Matsui 2019]
MATSUI, Takahiro, YAGASAKI, Takuma, MATSUMOTO, Masakazu and TANAKA, Hideki, 2019, Phase diagram of ice polymorphs under negative pressure considering the limits of mechanical stability. The Journal of Chemical Physics [online]. 28 January 2019. Vol. 150, no. 4, p. 041102. DOI 10.1063/1.5083021. Available from: http://dx.doi.org/10.1063/1.5083021

* [Matsumoto 2019]
MATSUMOTO, Masakazu, YAGASAKI, Takuma and TANAKA, Hideki, 2019, A Bayesian approach for identification of ice Ih, ice Ic, high density, and low density liquid water with a torsional order parameter. The Journal of Chemical Physics [online]. 7 June 2019. Vol. 150, no. 21, p. 214504. DOI 10.1063/1.5096556. Available from: http://dx.doi.org/10.1063/1.5096556

* [Matsumoto 2021]
MATSUMOTO, Masakazu, YAGASAKI, Takuma and TANAKA, Hideki, 2021, Novel Algorithm to Generate Hydrogen-Disordered Ice Structures. Journal of Chemical Information and Modeling [online]. 24 May 2021. DOI 10.1021/acs.jcim.1c00440. Available from: http://dx.doi.org/10.1021/acs.jcim.1c00440

* [Maynard-Casely 2010]
MAYNARD-CASELY, H. E., BULL, C. L., GUTHRIE, M., LOA, I., MCMAHON, M. I., GREGORYANZ, E., NELMES, R. J. and LOVEDAY, J. S., 2010, The distorted close-packed crystal structure of methane A. The Journal of Chemical Physics [online]. 14 August 2010. Vol. 133, no. 6, p. 064504. DOI 10.1063/1.3455889. Available from: http://dx.doi.org/10.1063/1.3455889

* [Mochizuki 2014]
MOCHIZUKI, Kenji, HIMOTO, Kazuhiro and MATSUMOTO, Masakazu, 2014, Diversity of transition pathways in the course of crystallization into ice VII. Phys. Chem. Chem. Phys. [online]. 2014. Vol. 16, no. 31, p. 16419–16425. DOI 10.1039/c4cp01616e. Available from: http://dx.doi.org/10.1039/c4cp01616e

* [Mousseau 2001]
MOUSSEAU, Normand and BARKEMA, G.T., 2001, Fast bond-transposition algorithms for generating covalent amorphous structures. Current Opinion in Solid State and Materials Science [online]. December 2001. Vol. 5, no. 6, p. 497–502. DOI 10.1016/s1359-0286(02)00005-0. Available from: http://dx.doi.org/10.1016/S1359-0286(02)00005-0

* [Nada 2003] Nada, Hiroki, and Jan P. J. M. van der Eerden. 2003. “An Intermolecular Potential Model for the Simulation of Ice and Water near the Melting Point: A Six-Site Model of H2O.” The Journal of Chemical Physics 118 (16): 7401–13.
* [Nakamura 2015]
NAKAMURA, Tatsuya, MATSUMOTO, Masakazu, YAGASAKI, Takuma and TANAKA, Hideki, 2015, Thermodynamic Stability of Ice II and Its Hydrogen-Disordered Counterpart: Role of Zero-Point Energy. The Journal of Physical Chemistry B [online]. 3 December 2015. Vol. 120, no. 8, p. 1843–1848. DOI 10.1021/acs.jpcb.5b09544. Available from: http://dx.doi.org/10.1021/acs.jpcb.5b09544

* [Petrenko 1999] Petrenko and Whitworth, Physics of Ice, Table 11.2.
* [Rosso 2016]
DEL ROSSO, Leonardo, CELLI, Milva and ULIVI, Lorenzo, 2016, New porous water ice metastable at atmospheric pressure obtained by emptying a hydrogen-filled ice. Nature Communications [online]. 7 November 2016. Vol. 7, no. 1. DOI 10.1038/ncomms13394. Available from: http://dx.doi.org/10.1038/ncomms13394

* [Russo 2014]
RUSSO, John, ROMANO, Flavio and TANAKA, Hajime, 2014, New metastable form of ice and its role in the homogeneous crystallization of water. Nature Materials [online]. 18 May 2014. Vol. 13, no. 7, p. 733–739. DOI 10.1038/nmat3977. Available from: http://dx.doi.org/10.1038/NMAT3977

* [Salzmann 2006]
SALZMANN, C. G., 2006, The Preparation and Structures of Hydrogen Ordered Phases of Ice. Science [online]. 24 March 2006. Vol. 311, no. 5768, p. 1758–1761. DOI 10.1126/science.1123896. Available from: http://dx.doi.org/10.1126/science.1123896

* [Sikiric 2010]
DUTOUR SIKIRIĆ, Mathieu, DELGADO-FRIEDRICHS, Olaf and DEZA, Michel, 2010, Space fullerenes: a computer search for new Frank–Kasper structures. Acta Crystallographica Section A Foundations of Crystallography [online]. 5 August 2010. Vol. 66, no. 5, p. 602–615. DOI 10.1107/s0108767310022932. Available from: http://dx.doi.org/10.1107/S0108767310022932

* [Smirnov 2013]
SMIRNOV, Grigory S. and STEGAILOV, Vladimir V., 2013, Toward Determination of the New Hydrogen Hydrate Clathrate Structures. The Journal of Physical Chemistry Letters [online]. 9 October 2013. Vol. 4, no. 21, p. 3560–3564. DOI 10.1021/jz401669d. Available from: http://dx.doi.org/10.1021/jz401669d

* [Stampfli 1986] Stampfli, P. A dodecagonal quasi-periodic lattice in 2 dimensions. Helv. Phys. Acta 59, 1260–1263 (1986).
* [Strobel 2016] Strobel, Timothy A et al. “Hydrogen-Stuffed, Quartz-Like Water Ice.” Journal of the American Chemical Society 138.42 (2016): 13786–13789.
* [Teeratchanan 2015]
TEERATCHANAN, Pattanasak and HERMANN, Andreas, 2015, Computational phase diagrams of noble gas hydrates under pressure. The Journal of Chemical Physics [online]. 21 October 2015. Vol. 143, no. 15, p. 154507. DOI 10.1063/1.4933371. Available from: http://dx.doi.org/10.1063/1.4933371

* [Vos 1993]
VOS, Willem L., FINGER, Larry W., HEMLEY, Russell J. and MAO, Ho-kwang, 1993, NovelH2-H2O clathrates at high pressures. Physical Review Letters [online]. 8 November 1993. Vol. 71, no. 19, p. 3150–3153. DOI 10.1103/physrevlett.71.3150. Available from: http://dx.doi.org/10.1103/PhysRevLett.71.3150

* [Weaire 1994]
WEAIRE, D. and FHELAN, R., 1994. The structure of monodisperse foam. Philosophical Magazine Letters [online]. November 1994. Vol. 70, no. 5, p. 345–350. DOI 10.1080/09500839408240997. Available from: http://dx.doi.org/10.1080/09500839408240997

* [Yagasaki 2018]
YAGASAKI, Takuma, MATSUMOTO, Masakazu and TANAKA, Hideki, 2018, Phase Diagrams of TIP4P/2005, SPC/E, and TIP5P Water at High Pressure. The Journal of Physical Chemistry B [online]. 17 July 2018. Vol. 122, no. 31, p. 7718–7725. DOI 10.1021/acs.jpcb.8b04441. Available from: http://dx.doi.org/10.1021/acs.jpcb.8b04441

* [Zhao 2019] Zhao, C.-L. et al. Seven-Site Effective Pair Potential for Simulating Liquid Water. J. Phys. Chem. B 123, 4594-4603 (2019).


# Algorithms and how to cite them.

The algorithms to make a depolarized hydrogen-disordered ice are explained in these papers:

M. Matsumoto, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered
Ice Generator",  J. Comput. Chem. 39, 61-64 (2017). [DOI: 10.1002/jcc.25077](http://doi.org/10.1002/jcc.25077)

    @article{Matsumoto:2017bk,
        author = {Matsumoto, Masakazu and Yagasaki, Takuma and Tanaka, Hideki},
        title = {GenIce: Hydrogen-Disordered Ice Generator},
        journal = {Journal of Computational Chemistry},
		volume = {39},
		pages = {61-64},
        year = {2017}
    }

M. Matsumoto, T. Yagasaki, and H. Tanaka, “Novel Algorithm to Generate Hydrogen-Disordered Ice Structures.”, J. Chem. Info. Modeling 61 (6): 2542–46 (2021). [DOI:10.1021/acs.jcim.1c00440](https://doi.org/10.1021/acs.jcim.1c00440)

    @article{Matsumoto:2021,
        author = {Matsumoto, Masakazu and Yagasaki, Takuma and Tanaka, Hideki},
        title = {Novel Algorithm to Generate Hydrogen-Disordered Ice Structures},
        journal = {Journal of Chemical Information and Modeling},
        volume = {61},
        pages = {2542-2546},
        year = {2021}
    }

# How to contribute

GenIce has been available as open source software on GitHub(https://github.com/vitroid/GenIce/) since 2015. Feedback, suggestions for improvements and enhancements, bug fixes, etc. are sincerely welcome. Developers and test users are also welcome. If you have any ice that is publicly available but not included in GenIce, please let me know.
