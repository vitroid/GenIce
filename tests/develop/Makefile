1.gro:
	../../genice.x 1h -s 1 > $@
2.gro:
	../../genice.x 1h -s 2 > $@
test.gro: 1.gro 2.gro
	cat 1.gro 1.gro 1.gro 1.gro 2.gro 2.gro 2.gro 2.gro > $@
test.yap: test.gro ../../genice/lattice.py
	../../analice.x test.gro -f yaplot > $@
