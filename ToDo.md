## Bug

* poetryはシンボリックリンクをパッケージに含めない。そのため、いくつかのモジュールが使えなくなる。
  * .genice2/を一時的に作成してsymlinkをコピーに変換し、そこからwheelを生成するようにした。

## ToDo

    
* Stage5が異常に遅い。なんで?
* hydrates with guests (using networkx.subgraph.isomorphism)
* semiclathrate hydrates
* hydrates of water and NH4F
* 各ステージが、独自の構造体を出力するようにする?
* 各ステージの入力と出力をもっと明確にする。全部引数で明示的に渡す。

## Future (GenIce3)
* Prefectによるワークフローの解体。ステージとかもう要らない。
* taipyかも。

## Done
* genice-core
* As a side effect, random module has been eliminated.
* Monoclinic cell for ice 5
* High pressure ices with prearranged network topology
* MDAnalysis integration.
* reshapeオプション自体はうまく動いた。ただ、a軸、b軸の向きをちゃんとしないとgromacs形式で書きだせない。
    * gromacsで書く場合には、grand_cellmatを早い段階で再調整する必要がある。
    * もとのセル自体も、gromacs仕様になっているとも限らない。
    * てっとり早いのは、a,b,c,A,B,C形式にして、戻すこと。
    * gromacs.pyだけの問題なので、そこに修正を加えた。