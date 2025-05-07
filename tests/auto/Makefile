# 改造したあとのbranchで実行する。
test:
	make -f Makefile.rules test

 
Makefile.rules: $(wildcard ../../genice2/*.py) mkmk.py
	python mkmk.py > $@

# 比較元(autotest-ref)のbranchで実行してファイルを準備する。
prepare: Makefile.rules
	make -f Makefile.rules -j 8 -k prepare 

clean:
	-rm *.gro Makefile.prepare


# * gro以外のフォーマット
# 以下はランダム選択でもいいかも
# * ion-doped
# * 乱数の違い
# * ゲスト、セミゲスト
# * プラグインのオプションごと。これについては、プラグインのほうでテスト方法を準備する必要がある。