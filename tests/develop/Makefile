all: 7.gro 7.re.gro 7.fail.gro
	echo Done.
7.gro:
	genice 7 -r 2 2 2 > $@
7.re.gro:
	genice gromacs[7.gro:O:H] -r 1 1 1 > $@
# It must fail because the molecular orientations are missing.
# ice 7 is so crowded structure that its HB network cannot be estimated solely from Oxygen positions.
7.fail.gro:
	genice gromacs[7.gro:O] -r 1 1 1 > $@

clean:
	-rm *.gro *~

