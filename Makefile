all:
	echo Hello.
%: temp_% genice/__main__.py
	genice -h | python3 Utilities/replace.py %%usage%% "    " $< > $@
%.rst: %.md
	md2rst $<

install:
	make README.rst
	./setup.py install
uninstall:
	-pip3 uninstall -y genice
pypi:
	make README.rst
	./setup.py check
	./setup.py sdist bdist_wheel upload
test:
	-mkdir trash
	./genice.x test1 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test2 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test3 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x test4 --guest 12=g12 --rep 1 1 1 --debug > /dev/null
	./genice.x 2d --rep 2 2 2 > trash/2d.gro
	./genice.x 3  --rep 2 2 2 --format c > trash/3.ar3a
	./genice.x 4  --rep 2 2 2 --format d > trash/4.ngph
	./genice.x 5  --rep 2 2 2 --format e > trash/5.nx3a
	./genice.x 6  --rep 2 2 2 --format g > trash/6.gro
	./genice.x 7  --rep 2 2 2 --format m > trash/7.mdv
	./genice.x 12 --rep 2 2 2 --format o > trash/12.scad
	./genice.x 16 --rep 1 1 1 --format p > trash/16.py
	./genice.x 17 --rep 1 1 1 --format q > trash/17.nx4a
	./genice.x 1c --rep 1 1 1 --format r > trash/1c.ar3r
	./genice.x 1h --rep 1 1 1 --format y > trash/1h.yap
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf GenIce.egg-info
	-rm README.rst
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf 
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
