CMD=lattices/_tool/convert_to_genice3.py


# shell scriptのけっか を変数に入れたい
ALL=$(shell cd lattices; ls *.py)

all: $(ALL)

%.py: lattices/%.py $(CMD)
	python $(CMD) $< .
