## ToDo

- [x] reactive なデータ生成。無駄がなくていい。
- [x] spot anion & cation
- 目玉機能。この 2 つが完成するまでは公表できない。
  - [ ] group!
  - [ ] topological defects
    - 指定方法を真剣に考える。本来なら、edgeを指定する必要がある。まあそれでいいか。-D 0 2 -L 4 5
    - GenIce3 CLIに登録する必要はないと思う。APIで簡単に書ければ十分。 
- [x] Loader (CIF, mdanalysis)
- [x] API レベルでの mdanalysis との連携方法。gro ファイルを経由するのが手っ取り早いのか。
- [ ] Exporter を増やす。
  - [x] svg 細かいチェックはできてない。
  - [x] png 細かいチェックはできてない。
  - [x] py3dmol
  - [x] mdanalysis これも py3dmol 方式で、一旦 gromacs を経由すれば簡単。コマンドラインで呼ぶ場合(dump)は、オプションには出力ファイル名を掻く。API で呼ぶ場合(universe())には、`Universe`オブジェクトを返す。
  - [x] \_KG まだおかしい。
- [x] 速度改善 →genice2 よりは早そう。
- [ ] Jupyter sample. xFAU はひとまず回避しよう。
- [ ] xFAU なぜうまくつながらないのか謎。
- [x] オプション内の Plugin の扱いの統一。コマンドラインではやむを得ないが、API から呼ぶ場合は loader を省略しないほうが良い。(オプションを指定する場合とで書式が変わってししまうのは結局親切といえない)
  - [x] 階層的オプションを扱うためのデータ形式の再検討。JSON5? HOCON?
  - [x] 深く再検討する。svg ですでにリスト型オプションのセパレータに`,`を使っていることに留意。
- [x] seed の設定は reactive か?
- [ ] Plugin を実行すると usage が表示されるように。
- [ ] help の充実。

## Bug

- poetry はシンボリックリンクをパッケージに含めない。そのため、いくつかのモジュールが使えなくなる。
  - .genice2/を一時的に作成して symlink をコピーに変換し、そこから wheel を生成するようにした。
