SEED=1234

$(SEED).makefile: make_makefile.py
	python make_makefile.py $(SEED) > $@


# 比較したいブランチを(手動で)チェックアウトしてから、make tag.[label]するとテスト実行してファイルがlabel/に作られる。
# これらを別のブランチでの結果と比較する。
tag.%: $(SEED).makefile
	# git checkout refs/tags/$*
	-mkdir $*
	make -f ../$(SEED).makefile -C $* -j 8 -k

compare: # tag.2.2.4 tag.2.2.5
	for file in `cd 2.2.4; ls`; do diff 2.2.4/$$file 2.2.5/$$file; done
