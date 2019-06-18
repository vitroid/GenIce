# Change log

## 1.0RC4

* GenIce no longer refers the files in the User global folder (.genice). Make them locally if necessary.
* The way to define cell dimension is changed.
* File loaders for analice tool are separated into a folder.

## 1.0RC3

* Several exotic ices are added.

## 1.0RC2

* Load and save multiple files (analice).

## 1.0RC1

* Documents are updated.

## 1.0RC

* Some functions for common use are separated into other packages. (PairList, etc)
* Some plugins that require special libraries are separated into extra packages (vpython, cif, zeolite, svg).
* New plugin hander that enables to implement extensions in separate packages.

## 0.24

* Random noise is no logger added to the molecular positions by default. (Use --add_noise option instead.)

## 0.22 (stable, release)

* Added AnalIce.

## 0.20.2, 0.20.3

* Added --version option.
* The version number is also shown in the usage.

## 0.20.1

* Atomic unit is supported in mdview format.

## 0.19 (develop), 0.20

* Added gromacs module as a lattice module in order to load a .gro
file as an ice structure.
* Added zeolite module as a lattice mofule.
* Added cif module.
* Added `--asis` option to use GenIce for file conversion.
* Changed the default lattice repetition numbers from [2,2,2] to [1,1,1]

## 0.18

* Direct graphical rendering with vpython.
* Added polygonnal expression in yaplot output.
* Added art examples for OpenSCAD format.

## 0.17

* svg_poly module.

## 0.16

* Ring phase statistics.
* Radial Kirkwood G function.
* Some plugins accept options using brackets.
* Cell reshaper.
* Added the current working path as a searchpath for the plugins.
* Aeroice generator.
* Accept the structure that does not obey the ice rules.

## 0.15.1

* Bug fix in case the atomic number exceeds 100 000 in Gromacs format.

## 0.15

* Simulation cell-related functions are separated to a module.

## 0.14

* Adapted to NetworkX 2.

## 0.13

* Regulated the value range of Euler's angles.

## 0.12

* Aeroices are added.

## 0.11

* Added "hooks" for the formats plugin.

## 0.10.6

* Added some new ice structures.

## 0.10.5

## 0.10.4

* Bug fix for ice 5.

## 0.10.3

* Bug fix in format modules.

## 0.10.2

* Bug fix.

## 0.10.1

* Accept semiclathrates, ion doping, etc.

## 0.10.0

* Accept hydrogen-ordered ices.

## 0.1

* First release. (Jun. 25, 2015)
