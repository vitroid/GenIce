GENICE=../../../genice.x
# .DELETE_ON_ERROR:


all: MTN.gro MTN2.gro
	echo Done.

# retrieve from IZA-structure DB
MTN.gro:
	$(GENICE) zeolite[MTN] -r 1 1 1 > $@

# Convert from a local CIF file
MTN2.gro: MTN.cif
	$(GENICE) cif[MTN.cif] -r 1 1 1 > $@

install:
	pip3 install cif2ice
clean:
	-rm *.gro *.cif
