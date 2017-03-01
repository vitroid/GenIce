all:
	echo Hello.
%: temp_% genice/__main__.py
	genice -h | python3 Utilities/replace.py %%usage%% "    " $< > $@
%.rst: %.md
	md2rst $<

install:
	make README.rst
	./setup.py install
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
test:
	genice 2d --rep 2 2 2 > 2d.gro
	genice 3  --rep 2 2 2 --format e > 3.nx3a
	genice 4  --rep 2 2 2 --format q > 4.nx4a
	genice 5  --rep 2 2 2 --format y > 5.yap
	genice 6  --rep 2 2 2 --format m > 6.mdv
	genice 7  --rep 2 2 2 --format d > 7.ngph
	genice 12 --rep 2 2 2 --format o > 12.scad
	genice 16 --rep 1 1 1 --format p > 16.py
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
