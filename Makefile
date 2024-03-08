PKGNAME=genice2
REPNAME=GenIce

all: README.md
	echo Hello.
%: temp_% Utilities/replacer.py genice2/__init__.py genice2/plugin.py citations.json pyproject.toml
	python Utilities/replacer.py < $< > $@

update-citations:
	cp citations.json old.citations.json
	python Utilities/citation.py < old.citations.json > citations.json
	-diff old.citations.json citations.json

test:
	make -C tests all

test-deploy:
	poetry publish --build -r testpypi
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PKGNAME)
uninstall:
	-pip uninstall -y $(PKGNAME)
build: README.md $(wildcard genice2/*.py)
	poetry build
deploy:
	poetry publish --build
check:
	poetry check




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
