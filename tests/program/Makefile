all:
	touch ref/*
	make -f Makefile.develop -j 8 -k
test_installed:
	make -f Makefile.installed -j 8
test_develop:
	make -f Makefile.develop

clean:
	-rm test* *~ @*
	-rm -rf sort_tmp lattices
	-rm *.cif *.test *.gro *.yap *.ar3a *.ringtest *.pov *.nx3a *.ar3a *.ar3r *.nx4a *.xyz *.mdv *.mdva *.py *.scad *.ngph *.rings *.kg *.ringstat *_gro *.png *.svg
