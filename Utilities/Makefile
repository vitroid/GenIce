C36.ar3a:
	~/github/TrajTools/xx 1 1 1 < ~/Dropbox/Tanaka/IceStructuresPro/2015/FKPhases/C36.ar3a > $@
	python3 ./ar3a2py.py < $@ > ~/Library/Application\ Support/GenIce/lattices/C36.py
	genice C36 --format y > C36.yap
	genice C36  > C36.gro
%.ar3a: ~/Dropbox/Tanaka/IceStructuresPro/2015/FKPhases/%.ar3a
	ln -s $< .
	python3 ./ar3a2py.py < $< > ~/Library/Application\ Support/GenIce/lattices/$*.py
	genice $* --format y > $*.yap
	genice $*  > $*.gro
#Cage positions from Dutour
%.ar3r: ~/Cluster/src/reatomiccoordinateofthespacefullerenes/MUfiles20/%.atom
	./sikiricparser.py < $< > $@

#1-decorated FAU
deco1.py:
	genice --rep 1 1 1 --format r 1c | ./decorator.py 1 > $@
test:
	./sikiricparser.py <~/Cluster/src/reatomiccoordinateofthespacefullerenes/MUfiles20/Struct13.atom
Struct%.yap:
	./sikiricparser.py <~/Cluster/src/reatomiccoordinateofthespacefullerenes/MUfiles20/Struct$*.atom > ~/Library/Application\ Support/GenIce/lattices/Struct$*.py
	genice Struct$* --rep 1 1 1 --format y > $@ 2> $@.log
all:
	i=1; while [ $$i -le 84 ]; do i=`printf '%02d' $$i`; echo Struct$$i.yap; i=`expr $$i + 1`; done | xargs make -j 8 -k 
clean:
	-rm deco1.py
