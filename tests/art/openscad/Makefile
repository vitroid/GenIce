GENICE=../../../genice.x


all:
	for s in ITT FAU RHO MTN sI; do echo $$s.scad; done | xargs make -k
clean:
	-rm *.scad *~

ITT.scad:
	$(GENICE) ITT -f openscad -r 3 3 2 > $@
FAU.scad:
	$(GENICE) FAU -f openscad -r 1 1 1 > $@
RHO.scad:
	$(GENICE) RHO -f openscad -r 1 1 1 > $@
MTN.scad:
	$(GENICE) 16 -f openscad -r 1 1 1 > $@
sI.scad:
	$(GENICE) sI -f openscad -r 1 1 1  > $@
