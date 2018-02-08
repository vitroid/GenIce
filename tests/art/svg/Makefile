GENICE=../../../genice.x


all:
	for s in ITT dtc FAU RHO MTN sI; do echo $$s.svg; done | xargs make -k
clean:
	-rm *.svg *~ formats/*~
	-rm -rf formats/__pycache__

ITT.svg: # formats/svg_poly.py Makefile
	$(GENICE) ITT -f svg_poly[-0.1,-0.1,-0.05] -r 3 3 2 > $@
dtc.svg: # formats/svg_poly.py Makefile
	$(GENICE) dtc -f svg_poly[-0.05,-0.05,-0.02]  -r 3 3 5 > $@
FAU.svg: # formats/svg_poly.py Makefile
	$(GENICE) FAU -f svg_poly -r 1 1 1 > $@
RHO.svg: # formats/svg_poly.py Makefile
	$(GENICE) RHO -f svg_poly -r 2 2 2 > $@
MTN.svg: # formats/svg_poly.py Makefile
	$(GENICE) 16 -f svg_poly[0.04,0.04,0.04] -r 2 2 2 > $@
sI.svg: # formats/svg_poly.py Makefile
	$(GENICE) sI -f svg_poly[-0.02,-0.02,-0.02] -r 2 2 2  > $@

prepare:
	pip3 install svgwrite
