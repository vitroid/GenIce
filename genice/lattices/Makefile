#MEP == CS1
#MTN == CS2, 16
#RHO == "sIII"
#BCT == "i"
CIF_FROM_ZEOLITEDB=MEP.cif MTN.cif RHO.cif BCT.cif
PY=$(patsubst %.cif, %.py, $(CIF_FROM_ZEOLITEDB))
all: $(CIF_FROM_EOLITEDB) $(PY)
%.py: %.cif
	./cifconv.py -i $< -o $@
BCT.py: BCT.cif
	./cifconv.py -i $< -o $@ -b 2 2 2
%.cif:
	curl http://asia.iza-structure.org/IZA-SC/cif/$@ > $@

prepare:
#	[ -e pycifrw ] || git clone https://bitbucket.org/jamesrhester/pycifrw.git
#	cd pycifrw && git checkout python3 && python3 setup.py install
#	brew install noweb
#I am sorry but the Python3 branch of pycifrw is not understandable for me (probably because the use of noweb)
#cifconv will work with python2 for a while.


