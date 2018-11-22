BASE=../../..

all: test1.yap test3.yap test5.yap

1.gro:
	$(BASE)/genice.x 1h -s 1 > $@
2.gro:
	$(BASE)/genice.x 1h -s 2 > $@
test.gro: 1.gro 2.gro
	cat 1.gro 1.gro 1.gro 1.gro 2.gro 2.gro 2.gro 2.gro > $@
test5.yap: test.gro $(BASE)/genice/lattice.py
	$(BASE)/analice.x test.gro --avgspan=5 -f yaplot > $@
test3.yap: test.gro $(BASE)/genice/lattice.py
	$(BASE)/analice.x test.gro --avgspan=3 -f yaplot > $@
test1.yap: test.gro $(BASE)/genice/lattice.py
	$(BASE)/analice.x test.gro --avgspan=1 -f yaplot > $@
clean:
	-rm *.yap *.gro
