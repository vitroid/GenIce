GENICE=../../genice.x
%.png: %.pov
	povray +I$< +W400 +H400 +D +FN +O$@ 

all: T.png
T.pov:
	$(GENICE) T -f povray -r 3 3 3 > T.pov
clean:
	rm T.*
