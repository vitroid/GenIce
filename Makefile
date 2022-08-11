all: README.md
	echo Hello.
%: temp_% Utilities/replacer.py genice2/__init__.py genice2/plugin.py citations.json
	-pip install jinja2
	python Utilities/replacer.py < $< > $@
	-fgrep '{{' $@

update-citations:
	cp citations.json old.citations.json
	python Utilities/citation.py < old.citations.json > citations.json
	-diff old.citations.json citations.json

pep8:
	autopep8 -r -a -a -i genice2/
test:
	make -C tests all
test-deploy: build
	-pip install twine
	twine upload -r pypitest dist/*
test-install: requirements.txt
	pip install -r $<
	pip install --index-url https://test.pypi.org/simple/ genice2


install:
	./setup.py install
uninstall:
	-pip uninstall -y genice2
build: $(wildcard genice2/*.py genice2/formats/*.py genice2/lattices/*.py genice2/molecules/*.py)
	./setup.py sdist # bdist_wheel


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
