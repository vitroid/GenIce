GENICE=../../genice.x
%.png: %.pov default.inc
	povray +I$< +W400 +H400 +D +FN +O$@ 

all: 1h.png
%.pov:
	$(GENICE) $* -f povray -r 3 3 3 > $@
clean:
	rm 1h.*
