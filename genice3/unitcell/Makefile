CMD=lattices/_tool/convert_to_genice3.py


# shell scriptのけっか を変数に入れたい
ALL=$(shell cd lattices; ls *.py)

all: $(ALL)

add: 
	for f in *.py; do grep genice2.lattices $$f || git add $$f; done

%.py: lattices/%.py $(CMD)
	python $(CMD) $< .
