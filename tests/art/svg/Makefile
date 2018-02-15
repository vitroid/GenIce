GENICE=../../../genice.x


all:
	for s in ITT FAU RHO MTN sI; do echo $$s.svg; done | xargs make -k
clean:
	-rm *.svg *~

ITT.svg:
	$(GENICE) ITT -f svg_poly[-0.1,-0.1,-0.05] -r 3 3 2 > $@
FAU.svg:
	$(GENICE) FAU -f svg_poly -r 1 1 1 > $@
RHO.svg:
	$(GENICE) RHO -f svg_poly -r 2 2 2 > $@
MTN.svg:
	$(GENICE) 16 -f svg_poly[0.04,0.04,0.04] -r 2 2 2 > $@
sI.svg:
	$(GENICE) sI -f svg_poly[-0.02,-0.02,-0.02] -r 2 2 2  > $@

prepare:
	pip3 install svgwrite
