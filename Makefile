all:
	echo Hello.
%: temp_% genice/__main__.py
	genice -h | python3 Utilities/replace.py %%usage%% "    " $< > $@
%.rst: %.md
	md2rst $<

install:
	make README.rst
	./setup.py install
uninstall:
	-pip3 uninstall -y genice
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
REFV=0.13
test:
	-mkdir trash
	./genice.x test1 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test2 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test3 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test4 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x 2d --rep 2 2 3 > trash/2d.gro                                   ;diff ref_$(REFV)/2d.gro trash/2d.gro
	./genice.x 3  --rep 2 2 3 --format c > trash/3.ar3a                        ;diff ref_$(REFV)/3.ar3a trash/3.ar3a
	./genice.x 4  --rep 2 2 3 --format d > trash/4.ngph                        ;diff ref_$(REFV)/4.ngph trash/4.ngph
	./genice.x 5  --rep 2 2 3 --format e > trash/5.nx3a                        ;diff ref_$(REFV)/5.nx3a trash/5.nx3a
	./genice.x 6  --rep 2 2 4 --format g > trash/6.gro                         ;diff ref_$(REFV)/6.gro  trash/6.gro
	./genice.x 7  --rep 2 2 4 --format m > trash/7.mdv                         ;diff ref_$(REFV)/7.mdv  trash/7.mdv
	./genice.x 12 --rep 2 2 3 --format o > trash/12.scad                       ;diff ref_$(REFV)/12.scad trash/12.scad
	./genice.x 16 --rep 1 1 2 --format p > trash/16.py                         ;diff ref_$(REFV)/16.py  trash/16.py
	./genice.x 17 --rep 1 1 2 --format q > trash/17.nx4a                       ;diff ref_$(REFV)/17.nx4a trash/17.nx4a
	./genice.x 1c --rep 1 1 2 --format r > trash/1c.ar3r                       ;diff ref_$(REFV)/1c.ar3r trash/1c.ar3r
#test for clathrates
	./genice.x CS1 --guest 14=g14*0.5 -G 0=me --rep 1 1 2 --format cif > trash/CS1.cif     ;diff ref_$(REFV)/CS1.cif trash/CS1.cif
#test for doped clathrate
	./genice.x CS2 --nodep -c 0=Na -a 1=Cl --rep 1 1 2 --format cif2 > trash/CS2.cif    ;Utilities/cifdiff.sh ref_$(REFV)/CS2.cif trash/CS2.cif
#test for spot semiclathrate (TBAB)
	./genice.x HS1 --nodep -c 3=N -a 1=Br -H 11=Bu-:3 -H 23=Bu-:3 -H 13=Bu-:3 -H 7=Bu-:3 --rep 1 1 2 --format xyz > trash/HS1.xyz; diff ref_$(REFV)/HS1.xyz trash/HS1.xyz
	./genice.x TS1 --guest 12=g12 --rep 1 1 2 --format yaplot > trash/TS1.yap  ;diff ref_$(REFV)/TS1.yap trash/TS1.yap
	#./genice.x T --guest 12=g12 --rep 1 1 1 --format povray > trash/T.pov  ;diff ref_$(REFV)/T.pov trash/T.pov
%.png: %.pov
	povray +I$< +W1000 +H1000 +D +FN +O$@ 
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
