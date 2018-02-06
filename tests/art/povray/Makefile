GENICE=../../../genice.x


all:
	for s in ITT dtc FAU RHO MTN; do echo $$s.pov $$s.png; done | xargs make -k
clean:
	-rm *.pov *.png *~ formats/*~
	-rm -rf formats/__pycache__

ITT.pov: # formats/povray_poly.py Makefile
	$(GENICE) ITT -f povray_poly[-0.1,-0.1,-0.05] -r 2 2 2 > $@
dtc.pov: # formats/povray_poly.py Makefile
	$(GENICE) dtc -f povray_poly[-0.05,-0.05,-0.02]  -r 3 3 5 > $@
FAU.pov: # formats/povray_poly.py Makefile
	$(GENICE) FAU -f povray_poly -r 1 1 1 > $@
RHO.pov: # formats/povray_poly.py Makefile
	$(GENICE) RHO -f povray_poly -r 2 2 2 > $@
MTN.pov: # formats/povray_poly.py Makefile
	$(GENICE) 16 -f povray_poly[0.04,0.04,0.04] -r 2 2 2 > $@

%.png: %.pov default.inc
	povray +I$< +W400 +H400 +D +FN +O$@ 
