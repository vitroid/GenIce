BASE=../..
GENICE=$(BASE)/genice.x


all:
	make -j 8 -k test_kg1.done test_kg2.done test_kg3.done test_kg4.done test_kg5.done test_kg6.done test_kg7.done

clean:
	-rm *.done *.log *~

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

test_kg%.done:
	i=1000; while [ $$i -lt 1100 ]; do echo kg$*_$$i.log; i=`expr $$i + 1`; done | xargs make -k
	cat kg$*_*.log | awk '{s[$$1]+=$$2;ss[$$1]+=$$2*$$2}END{for(x in s){a=s[x]/100;d=ss[x]/100-a*a;print x,a,sqrt(d/100)}}' | sort -n -k 1 > $@

