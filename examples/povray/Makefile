GENICE=../../genice.x
POVRAY=/opt/homebrew//Cellar/povray/3.7.0.10_2/bin/povray
LIBRARYPATH=/opt/homebrew//Cellar/povray/3.7.0.10_2/share/povray-3.7/include/
%.png: %.pov default.inc
	$(POVRAY) +L$(LIBRARYPATH) +I$< +W400 +H400 +D +FN +O$@ 

all: 1h.png
%.pov:
	$(GENICE) $* -f povray -r 3 3 3 > $@
clean:
	rm 1h.*
