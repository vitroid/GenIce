FILES=$(patsubst _original2/%, %, $(wildcard _original2/*.py))

all: $(FILES)

#%.py: _original/%.py _tool/__convert2.py
#	python _tool/__convert2.py $< > $@
%.py: _original2/%.py Makefile
	sed -e 's/genice /genice2 /g' -e 's/genice\./genice2./g' $< > $@
