# 改造したあとのbranchで実行する。
test:
	# run it in the new branch.
	python mkmk.py test > Makefile.test
	make -f Makefile.test

 
Makefile.prepare: $(wildcard ../../genice2/*.py)
	python mkmk.py > $@

# 比較元(autotest-ref)のbranchで実行してファイルを準備する。
prepare: Makefile.prepare
	make -f Makefile.prepare -j 8 -k 

clean:
	-rm *.gro


# * gro以外のフォーマット
# 以下はランダム選択でもいいかも
# * ion-doped
# * 乱数の違い
# * ゲスト、セミゲスト
# * プラグインのオプションごと。これについては、プラグインのほうでテスト方法を準備する必要がある。