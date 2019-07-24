all: README.md
	echo Hello.
%: temp_% genice/__main__.py genice/__init__.py
	./genice.x -h | python3 Utilities/replace.py %%usage%% $< > $@.1
	./analice.x -h | python3 Utilities/replace.py %%usage_analice%% $@.1 > $@




test:
	make -C tests all
test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install networkx numpy pairlist countrings yaplotlib
	pip install --index-url https://test.pypi.org/simple/ genice


install: README.md
	./setup.py install
uninstall:
	-pip uninstall -y genice
build: README.md $(wildcard genice/*.py genice/formats/*.py genice/lattices/*.py genice/molecules/*.py)
	./setup.py sdist bdist_wheel


deploy: build
	twine upload dist/*
check:
	./setup.py check
%.png: %.pov
	povray +I$< +W1000 +H1000 +D +FN +O$@ 
clean:
	-rm -rf build dist
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
