# 改造したあとのbranchで実行する。
test:
	make -f Makefile.rules test

 
Makefile.rules: $(wildcard ../../genice2/*.py) mkmk.py
	python mkmk.py > $@

# 比較元(autotest-ref)のbranchで実行してファイルを準備する。
prepare: Makefile.rules
	make -f Makefile.rules -k prepare 

clean:
	-rm *.output *.gro *.diff Makefile.rules
