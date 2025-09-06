![Logo]({{tool.genice.urls.logo}})

# GenIce2

{{tool.poetry.description}}

Version {{version}}

## New in GenIce2.2

- The core algorithm and its API are in a separate module, [`genice-core`](https://github.com/genice-dev/genice-core).

## Demo

The new GenIce works very well with interactive execution.
[Try instantly](https://colab.research.google.com/github/vitroid/GenIce/blob/main/jupyter.ipynb) on Google Colaboratory.

## Requirements

{% for item in tool.poetry.dependencies %}- {{item}}{{tool.poetry.dependencies[item]}}
{% endfor %}

<!-- **Note**: The package management system `poetry`, new in GenIce version 2.1, ignores all symlinks in package directories.
Because of this, some module "aliases" do not work correctly. (e.g. `genice2 1h` does not work, but `genice ice1h` does, because `1h.py` is an alias for `ice1h.py` .) -->

## Installation

GenIce is registered to [PyPI (Python Package Index)](https://pypi.python.org/pypi/GenIce).
Install with pip3.

```shell
% pip3 install genice2
```

## Uninstallation

```shell
% pip3 uninstall genice2
```

## Usage

{{usage}}

Use `./genice.x` instead of `genice2` if you want to use it inside the source tree.

## Examples

- To make a 3x3x3 units of a hydrogen-disordered ice IV (4) of TIP4P water in GROMACS
  .gro format:

```
          genice2 --water tip4p --rep 3 3 3  4 > ice4.gro
```

- To make a 2x2x4 units of CS2 clathrate hydrate structure of TIP4P water containing THF (united atom with a dummy site) in the large cage in GROMACS
  .gro format:

```
          genice2 -g 16=uathf6 --water tip4p --rep 2 2 4  CS2 > cs2-224.gro
```

## Basics

The program generates various ice lattice with proton disorder and without defect. The total dipole moment is always set to zero (except in the case you specify `--depol` option). The minimal structure (with --rep 1 1 1 option) is not always the unit cell of the lattice because it is difficult to deal with the hydrogen bond network topology of tiny lattice under periodic boundary condition. Note that the generated structure is not optimal according to the potential energy.

- To get a large repetition of ice Ih in XYZ format,

        genice2 --rep 8 8 8 1h --format xyz > 1hx888.xyz

- To get a ice V lattice of different hydrogen order in CIF format, use `-s` option to specify the random seed.

        genice2 5 -s 1024 --format cif > 5-1024.cif

- To obtain an ice VI lattice with different density and with TIP4P water model in gromacs format, use `--dens x` option to specify the density in g cm<sup>-3</sup>.

        genice2 6 --dens 1.00 --format g --water tip4p > 6d1.00.gro

GenIce is a modular program; it reads a unit cell data from a lattice plugin defined in the lattices folder, put water and guest molecules using a molecule plugin defined in the molecules/ folder, and output in various formats using a format plugin defined in the formats/ folder. You can write your plugins to extend GenIce. Some plugins also accept options.

## Clathrate hydrates

For clathrate hydrates, you can prepare the lattice with cages partially occupied by various guest molecules.

- To make a CS1 clathrate hydrate structure of TIP4P water containing CO2 in GROMACS .gro format: (60% of small cages are filled with co2 and 40% are methane)

        genice2 -g 12=co2*0.6+me*0.4 -g 14=co2 --water tip4p CS1 > cs1.gro

- To make a CS2 clathrate hydrate structure of TIP5P water containing THF molecules in the large cage, while only one cage is filled with methane molecule, first, just run `genice2` without guest specifications:

        genice2 CS2 > CS2.gro

  The list of cages will be output as follows:

        INFO   Cage types: ['12', '16']
        INFO   Cage type 12: {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183}
        INFO   Cage type 16: {136, 137, 138, 139, 140, 141, 142, 143, 16, 17, 18, 19, 20, 21, 22, 23, 160, 161, 162, 163, 164, 165, 166, 167, 40, 41, 42, 43, 44, 45, 46, 47, 184, 185, 186, 187, 188, 189, 190, 191, 64, 65, 66, 67, 68, 69, 70, 71, 88, 89, 90, 91, 92, 93, 94, 95, 112, 113, 114, 115, 116, 117, 118, 119}

  This indicates that there are two types of cages named `12` and `16`. Fill the `16` cages with THF and put a methane molecule in the `0`th cage of type `12` as follows:

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

        self.sites = np.array([[ 0.,    0., 0. ],
                               [-LCC/2, Y,  0. ],
                               [+LCC/2, Y,  0. ],])

        mass = np.array([16,14,14])
        # center of mass
        CoM = mass @ self.sites / np.sum(mass)
        self.sites -= CoM

        self.atoms_  = ["O","C","C"]
        self.labels_ = ["Oe","Ce","Ce"]
        self.name_   = "EO"
```

Write the code in eo.py. Make a folder named `molecules` in the current working directory and put it in.

_Note_: multiple occupancy is not implemented. If it is required, make a module of a virtual molecule that contains multiple molecules.

## Doping ions

Small ions may replace the host molecules. In that case, you can use `-a` and `-c` options to replace the specified water molecules with anions and cations.

The following example replaces the `0`th water molecule (in the replicated lattice) with Na cation and `1`st water molecule with Cl anion. The hydrogen bonds around the ions are organized appropriately.

    genice2 CS2 --depol=optimal -c 0=Na -a 1=Cl > CS2.gro

_Note 1_: The numbers of cations and anions must be the same. Otherwise, the ice rule is never satisfied and the program does not stop.

_Note 2_: The option `--depol=optimal` is also required because it is impossible to completely depolarize the structure containing ions.

_Note 3_: Protonic defects (H<sub>3</sub>O<sup>+</sup> and OH<sup>-</sup>) are not yet implemented.

## Semiclathrate hydrates

### Placement of a tetrabutylammonium ion (testing)

Let us assume that the id of the water molecule to be replaced by nitrogen of the TBA as zero. Place the nitrogen as a cation and also replace the water 2 by the counter-ion Br.

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

It indicates that the nitrogen is surrounded by cages with ids 0, 9, 2, and 7. Types for these cages can also be found in the info. Then, we put the Bu- group (minus does not mean ions) in these cages adjacent to dopant 0.

```shell
% genice2 HS1 -c 0=N -a 2=Br -H 0=Bu-:0 -H 9=Bu-:0 -H 2=Bu-:0 -H 7=Bu-:0 --depol=optimal > HS1.gro
```

Here the option `-H` specifies the group by `-H (cage id)=(group name):(root)`, and the root is the nitrogen that is specified by `-c` (cation) option.

### Placement of TBAB in the lattice module

_Under preparation_

It is more convenient if the lattice of the semiclathrate hydrate contains molecular ions in the appropriate locations in advance. Here we explain the way to make the special module for semclathrates.

## Output formats

| Name                     | Application                                                                                                           | extension    | water            | solute           | HB   | remarks                                                                                                                      |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------- | ------------ | ---------------- | ---------------- | ---- | ---------------------------------------------------------------------------------------------------------------------------- |
| `cif`, `cif2`            | CIF                                                                                                                   | `.cif`       | Atomic positions | Atomic positions | none | Experimental                                                                                                                 |
| `g`, `gromacs`           | [Gromacs](http://www.gromacs.org)                                                                                     | `.gro`       | Atomic positions | Atomic positions | none | Default format.                                                                                                              |
| `m`, `mdview`            | MDView                                                                                                                | `.mdv`       | Atomic positions | Atomic positions | auto |
| `mdv_au`                 | MDView                                                                                                                | `.mdv`       | Atomic positions | Atomic positions | auto | In atomic unit.                                                                                                              |
| `o`, `openscad`          | [OpenSCAD](http://www.openscad.org)                                                                                   | `.scad`      | Center of mass   | none             | o    | See tests/art/openscad for usage.                                                                                            |
| `povray`                 | Povray                                                                                                                | `.pov`       | Atomic positions | Atomic Positions | o    |
| `towhee`                 | TowHee                                                                                                                | `.coords`(?) | Atomic positions | Atomic positions | none |
| `xyz`                    | XYZ                                                                                                                   | `.xyz`       | Atomic positions | Atomic positions | none | Experimental                                                                                                                 |
| `exyz`                   | [extended XYZ](http://open-babel.readthedocs.io/en/latest/FileFormats/Extended_XYZ_cartesian_coordinates_format.html) | `.xyz`       | Atomic positions | Atomic positions | none | Extended XYZ format defined in Open Babel                                                                                    |
| `exyz2`                  | [extended XYZ](http://libatoms.github.io/QUIP/io.html#extendedxyz)                                                    | `.xyz`       | Atomic positions | Atomic positions | none | Extended XYZ format defined in QUIP                                                                                          |
| `y`, `yaplot`            | [Yaplot](https://github.com/vitroid/Yaplot)                                                                           | `.yap`       | Atomic positions | Atomic positions | o    | It renders molecular configurations and the HB network.                                                                      |
| `e`, `euler`             | Euler angles                                                                                                          | `.nx3a`      | Rigid rotor      | none             | none |
| `q`, `quaternion`        | Quaternions                                                                                                           | `.nx4a`      | Rigid rotor      | none             | none |
| `d`, `digraph`           | Digraph                                                                                                               | `.ngph`      | none             | none             | o    |
| `graph`                  | Graph                                                                                                                 | `.ngph`      | none             | none             | o    | Experimental.                                                                                                                |
| `c`, `com`               | CenterOfMass                                                                                                          | `.ar3a`      | Center of mass   | none             | none |
| `r`, `rcom`              | Relative CoM                                                                                                          | `.ar3r`      | Center of mass   | none             | none | In fractional coordinate system.                                                                                             |
| `p`, `python`, `reshape` | Python module                                                                                                         | `.py`        | Center of mass   | none             | none | Under development.                                                                                                           |
| `_ringstat`              | Ring phase statistics                                                                                                 |              |                  |                  |      | Statistical test suite 1: Check the appearance frequencies of the ring phases as a test for the intermediate-range disorder. |
| `rings`                  | [Yaplot](https://github.com/vitroid/Yaplot)                                                                           | `.yap`       | center of mass   | none             | o    | It renders HB rings.                                                                                                         |
| `_KG`                    | Kirkwood G(r)                                                                                                         |              |                  |                  |      | Statistical test suite 2: Calculate G(r) for checking long-range disorder in molecular orientations.                         |

By installing the [`genice2-mdanalysis`](https://github.com/vitroid/genice-mdanalysis) package separately, you can generate files in many formats for a large number of molecular dynamics package softwares. E.g.

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

| Symbol | <div style="width:300px">Description</div> |
| ------ | ------------------------------------------ |
{{ices}}

Ice names with double quotations are not experimentally verified.

You can prepare your own ice structures. Create a folder named `lattices` in the current working directory and put the plugins in it.

[cif2ice](https://github.com/vitroid/cif2ice) is a tool to retrieve a
cif file of zeolite from [IZA structure database](http://www.iza-structure.org/databases) and prepare a lattice
module in the path above.

Note: Different names are given in different nomenclature.

| CH/FI | CH  | ice | FK    | Zeo | Foam          |
| ----- | --- | --- | ----- | --- | ------------- |
| sI    | CS1 | -   | A15   | MEP | Weaire-Phelan |
| sII   | CS2 | 16  | C15   | MTN |               |
| sIII  | TS1 | -   | sigma | -   |               |
| sIV   | HS1 | -   | Z     | -   |               |
| sV    | HS2 | -   | \*    | -   |               |
| sVII  | CS4 | -   | \*    | SOD | Kelvin        |
| sH    | HS3 | -   | \*    | DOH |               |
| C0    | -   | 17  | \*    | -   |               |
| C1    | -   | 2   | \*    | -   |               |
| C2    | -   | 1c  | \*    | -   |               |

FI: Filled ices; CH: Clathrate hydrates; FK:Frank-Kasper duals; Zeo: Zeolites; Foam: foam crystals (Weaire 1994).

-: No correspondence; \*: Non-FK types.

Please ask [vitroid@gmail.com](mailto:vitroid@gmail.com) to add new ice structures.

## Water models

A water model can be chosen with `--water` option.

| symbol | type |
| ------ | ---- |
{{waters}}

## Guest molecules

| symbol | type |
| ------ | ---- |
{{guests}}

You can prepare your own guest molecules. Create a folder named `molecules` in the current working directory and put the plugins in it.

# Extra plugins

Some extra plugins are available via python package index using pip command.

For example, you can install RDF plugin by the following command,

```shell
% pip install genice2-rdf
```

And use it as an output format to get the radial distribution functions.

```shell
% genice2 TS1 -f _RDF > TS1.rdf.txt
```

## Output and analysis plugins

Analysis plugin is a kind of output plugin (specified with -f option).

| pip name                                                             | GenIce2 option         | Description                                                                                                         | output format                               | requirements              |
| -------------------------------------------------------------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------- | ------------------------------------------- | ------------------------- |
| [`genice2-cage`](https://github.com/vitroid/genice-cage)             | `-f _cage`             | Detect cages and quasi-polyhedra (vitrites).                                                                        | text, json, gromacs                         | `cycless`                 |
| [`genice2-rdf`](https://github.com/vitroid/genice-rdf)               | `-f _RDF`              | Radial distribution functions.                                                                                      | text                                        |                           |
| [`genice2-svg`](https://github.com/vitroid/genice-svg)               | `-f svg`<br />`-f png` | 2D graphics in SVG format.<br /> ... in PNG format.                                                                 | SVG<br />PNG                                | `svgwrite`                |
| [`genice2-twist`](https://github.com/vitroid/genice-twist)           | `-f twist`             | Calculate the twist order parameter (and visualize) [Matsumoto 2019]                                                | text (@BTWC)<br />SVG<br />PNG <br />yaplot | `twist-op`, `genice2-svg` |
| [`genice2-mdanalysis`](https://github.com/vitroid/genice-mdanalysis) | `-f mdanalysis`        | Output the atoms in various file formats that are served by [MDAnalysis](https://github.com/MDAnalysis/mdanalysis). | text, binary                                | `mdanalysis`              |

## Input plugins

Input plugins (a.k.a. lattice plugins) construct a crystal structure on demand.

| pip name                                               | GenIce2 usage                                       | Description                                                                       | requirements |
| ------------------------------------------------------ | --------------------------------------------------- | --------------------------------------------------------------------------------- | ------------ |
| [`genice2-cif`](https://github.com/vitroid/genice-cif) | `genice2 cif[ITT.cif]`<br /> `genice2 zeolite[ITT]` | Read a local CIF file as an ice structure.<br />Read a structure from Zeolite DB. | `cif2ice`    |

## New in GenIce2.1

GenIce2-MDAnalysis integration is now available. Try

```shell
% pip install genice2-mdanalysis
% genice2 1h -r 4 4 4 -f "mdanalysis[1h.pdb]"
```

to generate a PDB file.

## Changes from GenIce1

### Novel algorithm to make a structure obeying the ice rules in Stage 3.

- We have devised a completely new algorithm for orienting water molecules so that they obey ice rules. This algorithm can be applied only to defect-free ice. The algorithm runs in the following steps.
  1. First, based on the distances between neighboring molecules, the structure of the hydrogen-bond network is represented by a 4-connected undirected graph.
  2. The undirected graph is then randomly tiled with cycles. That is, we draw many cycles in the network so that all edges belong to only one of the cycles. It is always possible to a 4-connected regular graph.
  3. By directing each cycle, we can immediately obtain a directed graph that satisfies the ice rule. We can choose two orientations for each cycle so that the total polarization of the entire system is as small as possible.
  4. In rare cases, complete depolarization may not be possible. In such cases, it is depolarized in Stage 4.

### Faster, faster, faster.

Combinations of the new algorithm and other improvements in coding, the processing time of GenIce2 is about five times faster than that of GenIce1.

### Core algorithm is separated.

The core part of the new algorithm is separated as the TileCycles package.

### Colaboratory-ready!

Now GenIce2 works on the [Google Colaboratory!](https://colab.research.google.com/github/vitroid/GenIce/blob/genice2/jupyter.ipynb)

## New ices

Many new ice structures are added.

## Integration with MDAnalysis

GenIce2 is now integrated with MDAnalysis.

## References

{{citationlist}}

## Algorithms and how to cite them.

The algorithms to make a depolarized hydrogen-disordered ice are explained in these papers:

M. Matsumoto, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered
Ice Generator", J. Comput. Chem. 39, 61-64 (2017). [DOI: 10.1002/jcc.25077](http://doi.org/10.1002/jcc.25077)

    @article{Matsumoto:2017bk,
        author = {Matsumoto, Masakazu and Yagasaki, Takuma and Tanaka, Hideki},
        title = {GenIce: Hydrogen-Disordered Ice Generator},
        journal = {Journal of Computational Chemistry},
    	volume = {39},
    	pages = {61-64},
        year = {2017}
    }

M. Matsumoto, T. Yagasaki, and H. Tanaka, “GenIce-core: Efficient algorithm for generation of hydrogen-disordered ice structures.”, J. Chem. Phys. 160, 094101 (2024). [DOI:10.1063/5.0198056](https://doi.org/10.1063/5.0198056)

    @article{Matsumoto:2024,
        author = {Matsumoto, Masakazu and Yagasaki, Takuma and Tanaka, Hideki},
        title = {GenIce-core: Efficient algorithm for generation of hydrogen-disordered ice structures},
        journal = {Journal of Chemical Physics},
        volume = {160},
        pages = {094101},
        year = {2024}
    }

## How to contribute

GenIce has been available as open source software on GitHub({{project.urls.Homepage}}) since 2015.
Feedback, suggestions for improvements and enhancements, bug fixes, etc. are sincerely welcome.
Developers and test users are also welcome. If you have any ice that is publicly available but not included in GenIce, please let us know.
