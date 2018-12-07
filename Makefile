all: README.md
	echo Hello.
%: temp_% genice/__main__.py genice/__init__.py
	./genice.x -h | python3 Utilities/replace.py %%usage%% "    " $< > $@.1
	./analice.x -h | python3 Utilities/replace.py %%usage_analice%% "    " $@.1 > $@

test:
	make -C tests all
install:
	./setup.py install
uninstall:
	-pip uninstall -y genice
pypi:
	./setup.py check
	./setup.py sdist bdist_wheel upload
%.png: %.pov
	povray +I$< +W1000 +H1000 +D +FN +O$@ 
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
