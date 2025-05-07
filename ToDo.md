## Bug

- poetry はシンボリックリンクをパッケージに含めない。そのため、いくつかのモジュールが使えなくなる。
  - .genice2/を一時的に作成して symlink をコピーに変換し、そこから wheel を生成するようにした。

## ToDo

- GenIce3 に向けた整理を開始。
  - Hook を廃止し、Stage は内部でのみ利用するものとする。format plugin は自分が必要な要素を genice から取得する。genice は on demand でデータを生成する。
  - 改良状況は tests/auto/make で確認できる。
    - 現状、anion/cation ドーピングでひっかかっている。
    - ひととおり、現在の tests/auto/make を通るようになったら、commit する。でないと、中途半端すぎる。

* Stage5 が異常に遅い。なんで?
* hydrates with guests (using networkx.subgraph.isomorphism)
* semiclathrate hydrates
* hydrates of water and NH4F
* 各ステージが、独自の構造体を出力するようにする?
* 各ステージの入力と出力をもっと明確にする。全部引数で明示的に渡す。

## Future (GenIce3)

- Prefect によるワークフローの解体。ステージとかもう要らない。
- taipy かも。

## Done

- genice-core
- As a side effect, random module has been eliminated.
- Monoclinic cell for ice 5
- High pressure ices with prearranged network topology
- MDAnalysis integration.
- reshape オプション自体はうまく動いた。ただ、a 軸、b 軸の向きをちゃんとしないと gromacs 形式で書きだせない。
  - gromacs で書く場合には、grand_cellmat を早い段階で再調整する必要がある。
  - もとのセル自体も、gromacs 仕様になっているとも限らない。
  - てっとり早いのは、a,b,c,A,B,C 形式にして、戻すこと。
  - gromacs.py だけの問題なので、そこに修正を加えた。
