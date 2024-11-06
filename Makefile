PIPNAME=genice2
GITNAME=GenIce

all: README.md
	echo Hello.
%: temp_% Utilities/replacer.py genice2/__init__.py genice2/plugin.py citations.json pyproject.toml
	python Utilities/replacer.py < $< > $@


update-citations:
	cp citations.json old.citations.json
	python Utilities/citation.py < old.citations.json > citations.json
	-diff old.citations.json citations.json


# change symlinks to copies
prepare:
	-rm -rf .genice2
	mkdir .genice2
	rsync -avL --exclude="__pycache__" genice2/ .genice2/genice2/


test:
	make -C tests all

test-deploy: clean prepare build
	poetry publish -r testpypi
test-install:
	pip install --index-url https://test.pypi.org/simple/ $(PIPNAME)
uninstall:
	-pip uninstall -y $(PIPNAME)
build: README.md $(wildcard genice2/*.py) prepare
	poetry build -f wheel
deploy: clean prepare build
	poetry publish
check:
	poetry check


clone-myself-from-github:
	git clone -b genice-core https://github.com/vitroid/GenIce.git ../GenIce2.clone


%.png: %.pov
	povray +I$< +W1000 +H1000 +D +FN +O$@
clean:
	-rm -rf build dist .genice2
distclean:
	-rm *.scad *.yap @*
	-rm -rf build dist
	-rm -rf *.egg-info
	-rm .DS_Store
	find . -name __pycache__ | xargs rm -rf
	find . -name \*.pyc      | xargs rm -rf
	find . -name \*~         | xargs rm -rf
