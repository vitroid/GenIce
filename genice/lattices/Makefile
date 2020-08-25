FILES=$(patsubst _original/%, %, $(wildcard _original/*.py))
all: $(FILES)

%.py: _original/%.py _tool/__convert2.py
	python _tool/__convert2.py $< > $@
