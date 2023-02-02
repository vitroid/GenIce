GENICE=../../genice.x
validate: 1h.ringstat 2.ringstat 2d.ringstat 3.ringstat 4.ringstat 5.ringstat 6.ringstat 7.ringstat 9.ringstat 12.ringstat 16.ringstat
	echo Done.

# _ringstat module makes the statistics of bond orientations along each
# HB ring and compare the distribution with that of an ideal (isolated
# random) ring. The difference in the distribution is evaluated by
# Kullback-Leibler # divergence, d_{KL}.
# A typical dKL is zero for hydrogen-disordered ices, while it is
# larger than 1 for hydrogen-ordered ones like ice 2.
#
# Ringstat analysis validates the ring-scale randomness. GenIce tool
# also certifies the zero net dipole moment and Bernal-Fowler-Pauling
# ice rule in terms of the validity in global and local structures.

1h.%.ringstat:
	$(GENICE) 1h -r 24 24 24 -s $* -f _ringstat[max=6] > $@
%.ringstat:
	$(GENICE) $* -r 2 2 2  -f _ringstat[9] > $@
3.ringstat:
	$(GENICE) 3 -r 3 3 2  -f _ringstat[9] > $@
6.ringstat:
	$(GENICE) 6 -r 4 4 4  -f _ringstat[9] > $@
9.ringstat:
	$(GENICE) 9 -r 3 3 3  -f _ringstat[9] > $@





clean:
	-rm *~
	-rm *.ringstat
