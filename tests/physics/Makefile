BASE=../..
GENICE=$(BASE)/genice.x
REFD=ref

all: test_ring1.done test_ring2.done test_ring3.done test_kg1.done test_kg2.done test_kg3.done test_kg4.done test_kg5.done test_kg6.done test_kg7.done test_kg8.done test_kg9.done
	echo Done.

clean:
	-rm *.done *.log *~ */*~
	-rm -rf trash/



test_ring1.done:
	-mkdir trash
	$(GENICE) 1c --rep 8 8 8 --format _ring > trash/1c888.ringtest  ;diff $(REFD)/1c888.ringtest trash/1c888.ringtest
	touch $@
test_ring2.done:
	-mkdir trash
	$(GENICE) 2d --rep 4 2 8 --format _ring > trash/2d428.ringtest  ;diff $(REFD)/2d428.ringtest trash/2d428.ringtest
	touch $@
test_ring3.done:
	-mkdir trash
	$(GENICE) 3  --rep 4 4 4 --format _ring > trash/3x444.ringtest  ;diff $(REFD)/3x444.ringtest trash/3x444.ringtest
	touch $@



kg1_%.log:
	$(GENICE) 1c -s $* -f _KG -r 8 8 8  > $@
kg2_%.log:
	$(GENICE) 1c -s $* -f _KG -r 8 8 8 --nodep  > $@

kg3_%.log:
	$(GENICE) 2d -s $* -f _KG -r 4 2 8  > $@
kg4_%.log:
	$(GENICE) 2d -s $* -f _KG -r 4 2 8 --nodep  > $@
kg5_%.log:
	$(GENICE) 2  -s $* -f _KG -r 4 2 8  > $@

kg6_%.log:
	$(GENICE) 1h -s $* -f _KG -r 7 8 6  > $@
kg7_%.log:
	$(GENICE) 1h -s $* -f _KG -r 7 8 6 --nodep  > $@

# A test for radial Kirkwood G (KG).
# 1. Prepare a small (2x2x2) lattice ("semi-lattice") of hydrogen-disordered ice Ic.
# 2. Repeat it 4x4x4 times to make a larger lattice.
# 3. Calculate the averaged Kirkwood G for 100 different semi-lattices.

# Hydrogen-ordered semi-lattices
lattices/kg8_%.py:
	$(GENICE) 1c -s $* -f python -r 2 2 2 > $@
	$(GENICE) 1c -s $* -f fixed  -r 2 2 2 >> $@  #Fix the HB orientations.
                                                    #"fixed" is a local plugin just for this test.
# _KG calculate the radial KG function
kg8_%.log: lattices/kg8_%.py
	$(GENICE) kg8_$* -f _KG -r 4 4 4 > $@



lattices/kg9_%.py:
	$(GENICE) 1c -s $* -f python -r 4 4 4 > $@
	$(GENICE) 1c -s $* -f fixed  -r 4 4 4 >> $@  #Fix the HB orientations.
                                                    #"fixed" is a local plugin just for this test.
# Repeat the semi-lattices.
kg9_%.log: lattices/kg9_%.py
	$(GENICE) kg9_$* -f _KG -r 2 2 2 > $@

test_kg%.done:
	-mkdir trash
	i=1000; while [ $$i -lt 1100 ]; do echo kg$*_$$i.log; i=`expr $$i + 1`; done | xargs make -k
	cat kg$*_*.log | awk '{s[$$1]+=$$2;ss[$$1]+=$$2*$$2}END{for(x in s){a=s[x]/100;d=ss[x]/100-a*a;print x,a,sqrt(d/100)}}' | sort -n -k 1 > trash/$@
	diff $(REFD)/$@ trash/$@
	rm kg$*_*.log;	ln -s trash/$@ .


# Prepare a small cell of hydrogen-disordered.
lattices/kg10_%.py:
	$(GENICE) 1c -s $* -f python -r 2 2 2 > $@
	$(GENICE) 1c -s $* -f fixed  -r 2 2 2 >> $@  #Specify the HB orientations.

# Tile them and shuffle it by cycle inversions.
kg10_%.log: lattices/kg10_%.py
	$(GENICE) kg10_$* -f inversion -r 4 4 4 > $@
