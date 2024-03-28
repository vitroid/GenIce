GENICE=../../genice.x


all: prepare
	pytest test_format.py
	pytest test_lattice.py
	pytest test_loader.py
	pytest test_molecule.py

prepare: 3.ar3r 5.exyz 7.gro CS1.mdv CS2.mdva HS1.nx3a

%.ar3r: Makefile
	$(GENICE) $* -w tip4p -f rcom > $@
%.exyz: Makefile
	$(GENICE) $* -w tip4p -f exyz > $@
%.gro: Makefile
	$(GENICE) $* -w tip4p -f gromacs > $@
%.mdv: Makefile
	$(GENICE) $* -w tip4p -f mdview > $@
%.mdva: Makefile
	$(GENICE) $* -w tip4p -f mdv_au > $@
%.nx3a: Makefile
	$(GENICE) $* -w tip4p -f euler > $@


