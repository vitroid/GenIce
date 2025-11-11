## Bug

## ToDo

### Topological defect

- 一応、cation のみ、anion のみでも動くようにする。
- D,L の場合は、あとで水素の付加/除去を行って H3O と OH を欠陥に変換する。
- 問題はそれをコマンドラインでどうやって指示するか。
- というか、今後のことを考えるとコマンドラインではなく関数としてどう記述するかをまず考えよ。

  - eigenize()と hydroxify()を作る。事前に辺の向きを固定する。
    - Stage X に追加処理を加える、という方法がない。formats プラグインにするのは気持ち悪い。
    - 従来通りにやるなら、GenIce クラスの生成時のオプションを増やすしかない。
    - modular program と言いながら、ごっつい monolith だな。GenIce3 では解決したい。
      - format も追加機能も全部単一形状の plugin class を継承する。
      - 旧 panojector のように、plugin 自身がオプションを受けつける。
  - グラフの生成のあと、分子を置く段になって、特殊処理を行う。
    - Bjerrum の場合も同じ場所で分子を置く。

- Stage5 が異常に遅い。なんで?
- hydrates with guests (using networkx.subgraph.isomorphism)
- semiclathrate hydrates
- hydrates of water and NH4F
- 各ステージが、独自の構造体を出力するようにする?
- 各ステージの入力と出力をもっと明確にする。全部引数で明示的に渡す。

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
