# genice2.lattices → genice3.unitcell 変換スクリプト

## 概要

`convert_to_genice3.py`は、genice2 の lattices モジュールを genice3 の unitcell モジュールに自動変換するスクリプトです。

## 使い方

### 単一ファイルの変換

```bash
python3 genice2/lattices/_tool/convert_to_genice3.py <lattice_file.py> [output_dir]
```

例:

```bash
python3 genice2/lattices/_tool/convert_to_genice3.py genice2/lattices/RHO.py
# -> genice3/unitcell/RHO.py に出力

python3 genice2/lattices/_tool/convert_to_genice3.py genice2/lattices/A15.py /tmp/output
# -> /tmp/output/A15.py に出力
```

### すべてのファイルを一括変換

```bash
python3 genice2/lattices/_tool/convert_to_genice3.py --all [output_dir]
```

例:

```bash
python3 genice2/lattices/_tool/convert_to_genice3.py --all
# -> genice3/unitcell/ にすべて出力

python3 genice2/lattices/_tool/convert_to_genice3.py --all /tmp/output
# -> /tmp/output/ にすべて出力
```

## 対応パターン

### 1. 標準パターン（pairs 文字列あり）

```python
# genice2/lattices/A15.py
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.pairs = """
        0 1
        1 2
        """
        self.waters = """
        0.5 0.5 0.5
        """
        self.cell = cellvectors(a=1.0, b=1.0, c=1.0)
```

→ `pairs`を`nx.Graph`に変換して`graph`パラメータに渡します。

### 2. pairs なしパターン（自動生成）

```python
# genice2/lattices/RHO.py
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.waters = """
        0.5 0.5 0.5
        """
        self.bondlen = 3
        self.cell = cellvectors(a=1.0, b=1.0, c=1.0)
```

→ `graph=None`で自動生成されます。

### 3. CIF パターン

```python
# genice2/lattices/ice2.py
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.waters, self.fixed = CIF.waters_and_pairs(...)
        self.pairs = self.fixed
```

→ **警告を出して手動変換を促します**。CIF パターンは複雑なため、自動変換は不完全です。

### 4. fixed パターン

```python
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.fixed = """
        0 1
        """
```

→ `fixed`を`nx.DiGraph`に変換して`fixed`パラメータに渡します。

## 変換結果の確認

変換後は、以下の点を確認してください：

1. **CIF パターン**: 警告が出た場合は手動で変換が必要です
2. **cell 定義**: 複雑な計算式がある場合は手動で確認が必要な場合があります
3. **density**: コメントアウトされていますが、必要に応じて有効化してください
4. **fixed**: `self.pairs = self.fixed`のようなパターンは適切に変換されます

## 注意事項

- 変換スクリプトは基本的なパターンに対応していますが、すべてのケースをカバーするわけではありません
- CIF パターンや複雑な計算式を含む場合は、手動での確認・修正が必要です
- 変換後は必ずテストを実行して、正しく動作することを確認してください

## 例

### RHO.py の変換例

**変換前** (genice2/lattices/RHO.py):

```python
class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.waters = """..."""
        self.coord = "absolute"
        self.bondlen = 3
        self.density = 0.604398971981
        self.cell = cellvectors(a=13.339813507, ...)
```

**変換後** (genice3/unitcell/RHO.py):

```python
class UnitCell(genice3.unitcell.UnitCell):
    def __init__(self, **kwargs):
        graph = None  # pairsがない場合は自動生成
        waters = np.fromstring("""...""", sep=" ").reshape(-1, 3)
        coord = "absolute"
        bondlen = 3
        # density = 0.604398971981
        cell = cellvectors(a=13.339813507, ...)
        super().__init__(
            cell=cell,
            waters=waters,
            coord=coord,
            bondlen=bondlen,
            # density=density,
            **kwargs,
        )
```
