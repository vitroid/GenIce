all: 3.compare 4.compare
%.compare: %flat.sh.out %.py.out %.yaml.out %.sh.out
	-diff $*.py.out $*.yaml.out
	-diff $*.py.out $*flat.sh.out
	-diff $*.py.out $*.sh.out
	-diff $*.py.log $*flat.sh.log
	-diff $*.py.log $*.yaml.log
	-diff $*.py.log $*.sh.log


%.sh.out: %.sh Makefile
	sh $< > $@ 2> $<.log
%.py.out: %.py Makefile
	python3 $< > $@ 2> $<.log
%.yaml.out: %.yaml Makefile
	python3 -m genice3.cli.genice --config $< > $@ 2> $<.log

clean:
	\rm *.out
