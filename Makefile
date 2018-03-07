all:
	echo Hello.
%: temp_% genice/__main__.py genice/__init__.py
	./genice.x -h | python3 Utilities/replace.py %%usage%% "    " $< > $@
%.rst: %.md
	md2rst $<
test:
	make -C tests all
install:
	make README.rst
	./setup.py install
uninstall:
	-pip3 uninstall -y genice
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
%.png: %.pov
	povray +I$< +W1000 +H1000 +D +FN +O$@ 
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
