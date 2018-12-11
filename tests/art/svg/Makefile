GENICE=../../../genice.x


all:
	for s in ITT FAU RHO MTN sI; do echo $$s.svg; done | xargs make -k
clean:
	-rm *.svg *~

ITT.svg:
	$(GENICE) ITT -f svg[shadow] -r 3 3 2 > $@
FAU.svg:
	$(GENICE) FAU -f svg[shadow] -r 1 1 1 > $@
RHO.svg:
	$(GENICE) RHO -f svg[shadow] -r 2 2 2 > $@
MTN.svg:
	$(GENICE) 16 -f svg[shadow] -r 2 2 2 > $@
sI.svg:
	$(GENICE) sI -f svg[shadow] -r 2 2 2  > $@

prepare:
	pip3 install svgwrite
