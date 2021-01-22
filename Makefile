all: README.md
	echo Hello.
%: temp_% Utilities/replacer.py genice2/__init__.py genice2/plugin.py
	python Utilities/replacer.py < $< > $@
	-fgrep '{{' $@

update-citations:
	cp citations.json old.citations.json
	python Utilities/citation.py < old.citations.json > citations.json


test:
	make -C tests all
test-deploy: build
	twine upload -r pypitest dist/*
test-install:
	pip install networkx numpy pairlist countrings yaplotlib
	pip install --index-url https://test.pypi.org/simple/ genice2


install: README.md
	./setup.py install
uninstall:
	-pip uninstall -y genice2
build: README.md $(wildcard genice2/*.py genice2/formats/*.py genice2/lattices/*.py genice2/molecules/*.py)
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
	-rm -rf *.egg-info
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
