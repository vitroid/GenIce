BASE=../..
GENICE=$(BASE)/genice.x
LIB=$(BASE)/genice
FMT=$(LIB)/formats
UTIL=.
REFD=ref
.DELETE_ON_ERROR:

all: test1.done   test2.done   test3.done   test4.done   test_2d.gro
all: test_3.ar3a  test_4.ngph  test_5.nx3a  test_6.gro   test_7.mdv
all: test_12.scad test_16.py   test_17.nx4a test_1c.ar3r test_CS1.cif
all: test_CS2.cif test_HS1.xyz test_TS1.yap test_2dL.gro test_7.mdva
#all: test19.done
all: test_1hreshape.py test_2D3.ringtest test_prism7.ringtest
all: test_3xFAU.yap
all: test_T.pov
all: test_7.re.gro
all: test_MTN.gro test_MTN.re.gro
all: test_ITT.gro

clean:
	-rm test* *~
	-rm -rf sort_tmp
	-rm *.cif


test1.done: $(FMT)/gromacs.py
	$(GENICE) test1 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	touch test1.done
test2.done: $(FMT)/gromacs.py
	$(GENICE) test2 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	touch test2.done
test3.done: $(FMT)/gromacs.py
	$(GENICE) test3 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	touch test3.done
test4.done: $(FMT)/gromacs.py
	$(GENICE) test4 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	touch test4.done
test_2d.gro: $(FMT)/gromacs.py
	$(GENICE) 2d --rep 2 2 3 > $@                                   ;diff $(REFD)/$@ $@
test_3.ar3a: $(FMT)/com.py
	$(GENICE) 3  --rep 2 2 3 --format c > $@                        ;diff $(REFD)/$@ $@
test_4.ngph: $(FMT)/digraph.py
	$(GENICE) 4  --rep 2 2 3 --format d > $@                        ;diff $(REFD)/$@ $@
test_5.nx3a:  $(FMT)/euler.py
	$(GENICE) 5  --rep 2 2 3 --format e > $@                        ;diff $(REFD)/$@ $@
test_6.gro: $(FMT)/gromacs.py
	$(GENICE) 6  --rep 2 2 4 --format g > $@                        ;diff $(REFD)/$@ $@
test_7.mdv: $(FMT)/mdview.py
	$(GENICE) 7  --rep 2 2 4 --format m > $@                        ;diff $(REFD)/$@ $@
test_7.mdva: $(FMT)/mdv_au.py
	$(GENICE) 7  --rep 2 2 4 --format mdv_au > $@                   ;diff $(REFD)/$@ $@
test_12.scad: $(FMT)/openscad.py
	$(GENICE) 12 --rep 2 2 3 --format o > $@                        ;diff $(REFD)/$@ $@
test_16.py: $(FMT)/python.py
	$(GENICE) 16 --rep 1 1 2 --format p > $@                        ;diff $(REFD)/$@ $@
test_17.nx4a: $(FMT)/quaternion.py
	$(GENICE) 17 --rep 1 1 2 --format q > $@                        ;diff $(REFD)/$@ $@
test_1c.ar3r: $(FMT)/rcom.py
	$(GENICE) 1c --rep 1 1 2 --format r > $@                        ;diff $(REFD)/$@ $@
test_CS1.cif: $(FMT)/cif.py
	$(GENICE) CS1 --guest 14=g14*0.5 -G 0=me --rep 1 1 2 --format cif > $@ ;diff $(REFD)/$@ $@
test_CS2.cif: $(FMT)/cif2.py
	$(GENICE) CS2 --nodep -c 0=Na -a 1=Cl --rep 1 1 2 --format cif2 > $@  ;$(UTIL)/cifdiff.sh $(REFD)/$@ $@
test_HS1.xyz: $(FMT)/xyz.py
#test for spot semiclathrate (TBAB)
	$(GENICE) HS1 --nodep -c 3=N -a 1=Br -H 11=Bu-:3 -H 23=Bu-:3 -H 13=Bu-:3 -H 7=Bu-:3 --rep 1 1 2 --format xyz > $@; diff $(REFD)/$@ $@
test_TS1.yap: $(FMT)/yaplot.py
	$(GENICE) TS1 --guest 12=g12 --rep 1 1 2 --format yaplot > $@ ;diff $(REFD)/$@ $@
test_T.po: $(FMT)/povray.py
	$(GENICE) T --guest 12=g12 --rep 1 1 1 --format povray > $@   ;diff $(REFD)/$@ $@
test_2dL.gro: $(FMT)/gromacs.py
	$(GENICE) 2d --rep 8 8 8 > $@
test_1hreshape.py: $(FMT)/reshape.py
	$(GENICE) 1h --rep 1 1 1 --format reshape[1,1,-2,1,-1,0,1,1,1] > $@ ;diff $(REFD)/$@ $@
test_2D3.ringtest: $(FMT)/_ring.py
	$(GENICE) 2D3 --rep 5 5 1 --format _ring > $@                ;diff $(REFD)/$@ $@
test_prism7.ringtest: $(FMT)/_ring.py
	$(GENICE) prism[7] --rep 1 1 1 --format _ring > $@           ;diff $(REFD)/$@ $@
test_3xFAU.yap: $(FMT)/yaplot.py
	$(GENICE) xFAU[3] --rep 1 1 1 --format yaplot > $@           ;diff $(REFD)/$@ $@
test_7.gro: $(FMT)/gromacs.py
	$(GENICE) 7 -r 2 2 2 > $@
test_7.re.gro: test_7.gro $(FMT)/gromacs.py
	$(GENICE) gromacs[$<:O:H] --rep 1 1 1 --asis  > $@           ;diff $(REFD)/$@ $@
	-$(GENICE) gromacs[$<:O]  --rep 1 1 1 > $@ # should fail
	echo It must fail because the molecular orientations are missing.  Ice 7 is so crowded structure that its HB network cannot be estimated solely from Oxygen positions.
test_MTN.gro: $(FMT)/gromacs.py
	$(GENICE) zeolite[MTN] --rep 1 1 1 > $@
test_MTN.re.gro: MTN.cif $(FMT)/gromacs.py
	$(GENICE) cif[MTN.cif] --rep 1 1 1 > $@                     ;diff $(REFD)/$@ $@
test_ITT.gro: $(FMT)/gromacs.py
	$(GENICE) zeolite[ITT] --rep 1 1 1 > $@                     ;diff $(REFD)/$@ $@
