SEED=1234

$(SEED).makefile: make_tests.py
	python make_tests.py $(SEED) > $@

tag.%: $(SEED).makefile
	git checkout refs/tags/$*
	-mkdir $*
	make -f ../$(SEED).makefile -C $*

compare: # tag.2.2.4 tag.2.2.5
	for file in `cd 2.2.4; ls`; do diff 2.2.4/$$file 2.2.5/$$file; done
