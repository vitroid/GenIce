GENICE=../../genice.x
# validate: 1h.ringstat 2.ringstat 2d.ringstat 3.ringstat 4.ringstat 5.ringstat 6.ringstat 7.ringstat 9.ringstat 12.ringstat 16.ringstat
	# echo Done.

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

test: test0 test1 test2
test0:
	for n in 2 4; do \
	    rep=1000; \
		for i in `seq 1 $$rep`; do echo stats/$${n}_$$i.ringstat; done; \
	done | xargs make -j 8 -k
test1:
	for n in 8 16; do \
	    rep=100; \
		for i in `seq 1 $$rep`; do echo stats/$${n}_$$i.ringstat; done; \
	done | xargs make -j 8 -k
test2:
	for n in 32; do \
	    rep=32; \
		for i in `seq 1 $$rep`; do echo stats/$${n}_$$i.ringstat; done; \
	done | xargs make 


%.ringstat:
	n=`basename $* | sed -e 's/_.*//'` ; \
	s=`basename $* | sed -e 's/^.*_//'`; \
	echo $$n $$s; \
	$(GENICE) 1c -r $$n $$n $$n -f _ringstat -s $$s > $@



# 1h.%.ringstat:
# 	$(GENICE) 1h -r 24 24 24 -s $* -f _ringstat[max=6] > $@
# %.ringstat:
# 	$(GENICE) $* -r 2 2 2  -f _ringstat[9] > $@
# 3.ringstat:
# 	$(GENICE) 3 -r 3 3 2  -f _ringstat[9] > $@
# 6.ringstat:
# 	$(GENICE) 6 -r 4 4 4  -f _ringstat[9] > $@
# 9.ringstat:
# 	$(GENICE) 9 -r 3 3 3  -f _ringstat[9] > $@





clean:
	-rm *~
	-rm *.ringstat
