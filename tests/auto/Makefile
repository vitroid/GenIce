test:
	# run it in the new branch.
	python mkmk.py test > Makefile.test
	make -f Makefile.test

Makefile.prepare:
	python mkmk.py > $@

prepare: Makefile.prepare
	# run it in the autotest-ref branch.
	make -f Makefile.prepare
