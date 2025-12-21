3.compare: 3flat.sh.out 3.py.out 3.yaml.out 3.sh.out
	-diff 3.py.out 3.yaml.out
	-diff 3.py.out 3flat.sh.out
	-diff 3.py.out 3.sh.out


%.sh.out: %.sh Makefile
	sh $< > $@ 2> $<.log
%.py.out: %.py Makefile
	python3 $< > $@ 2> $<.log
%.yaml.out: %.yaml Makefile
	python3 -m genice3.cli.genice --config $< > $@ 2> $<.log

clean:
	\rm *.out
