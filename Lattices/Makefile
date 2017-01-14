CIF_FROM_ZEOLITEDB=MEP.cif MTN.cif RHO.cif
PY=$(patsubst %.cif, %.py, $(CIF_FROM_ZEOLITEDB))
all: $(CIF_FROM_EOLITEDB) $(PY)
%.py: %.cif
	./cifconv.py -i $< -o $@
%.cif:
	curl http://asia.iza-structure.org/IZA-SC/cif/$@ > $@

