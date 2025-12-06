# GenIce3 の Dagster への移行設計

## 概要

現在の`GenIce3`クラスは、`DependencyCacheMixin`と`property_depending_on`デコレータを使用して、reactive な依存関係管理を実装しています。これを Dagster の asset システムに移行することで、より標準的で保守しやすいコードベースにすることができます。

## 設計方針

### 1. クラスの解体

**結論: GenIce3 クラスは解体すべきです。**

理由:

- **関数型プログラミングの原則**: クラス変数のような内部状態は、関数型プログラミングでは「諸悪の根源」とみなされます。状態を持つことで、テストが困難になり、並列実行が難しくなります。
- **Dagster の設計思想**: Dagster は asset ベースのデータパイプラインを推奨しており、各 asset は独立した関数として実装されます。
- **明示的な依存関係**: クラス内の暗黙的な依存関係ではなく、asset 依存関係として明示的に宣言することで、コードの可読性と保守性が向上します。

### 2. 設定の管理

設定パラメータ（`depol_loop`, `replication_matrix`, `unitcell`）は以下の方法で管理できます:

#### 方法 A: @dataclass を使用（推奨）

```python
@dataclass
class GenIce3Config:
    depol_loop: int = 1000
    replication_matrix: np.ndarray = None
    unitcell: UnitCell = None
```

#### 方法 B: Dagster の Config を使用

```python
class GenIce3Config(Config):
    depol_loop: int = Field(default=1000, description="Depolarization loop")
    replication_matrix: Optional[List[List[int]]] = Field(default=None)
```

#### 方法 C: AssetKey を使用

設定を asset として定義し、他の asset が依存するようにする。

### 3. 計算結果の asset 化

現在の`@property_depending_on`でデコレートされたプロパティは、すべて独立した asset 関数に変換します:

```python
# 現在の実装
@property_depending_on("replication_matrix", "unitcell")
def cell(self):
    return self.unitcell.cell @ self.replication_matrix

# Dagster版
@asset
def cell(config: GenIce3Config, unitcell: UnitCell) -> np.ndarray:
    return unitcell.cell @ config.replication_matrix
```

### 4. 依存関係の宣言

Dagster では、asset 間の依存関係を明示的に宣言します:

```python
@asset
def graph(
    unitcell: UnitCell,
    replica_vectors: np.ndarray,
    replica_vector_labels: Dict[Tuple[int, int, int], int],
    config: GenIce3Config,
) -> nx.Graph:
    # ...
```

## 移行のメリット

1. **テスト容易性**: 各 asset は独立した関数なので、単体テストが容易です。
2. **並列実行**: Dagster が自動的に依存関係を解析し、並列実行可能な asset を並列実行します。
3. **可視化**: Dagster UI で依存関係グラフを可視化できます。
4. **再現性**: 各 asset の実行結果をキャッシュし、再実行を避けることができます。
5. **デバッグ**: 特定の asset だけを再実行することができます。

## 移行の課題

1. **既存コードとの互換性**: 既存のコードが`GenIce3`クラスを使用している場合、移行が必要です。
2. **設定の管理**: 設定をどのように管理するか（Config、AssetKey、dataclass）を決定する必要があります。
3. **メソッドの変換**: `water_molecules`、`guest_molecules`などのメソッドも関数に変換する必要があります。

## 実装例

詳細な実装例は`genice2/genice3_dagster.py`を参照してください。

## 移行手順

1. **段階的移行**: まず、新しい asset ベースの実装を作成し、既存のクラスベースの実装と並行して動作させる。
2. **テスト**: 両方の実装で同じ結果が得られることを確認する。
3. **段階的置き換え**: 既存のコードを段階的に新しい実装に置き換える。
4. **削除**: すべてのコードが新しい実装を使用するようになったら、古いクラスベースの実装を削除する。

## 結論

GenIce3 クラスは解体し、すべてを関数と asset に変換することを推奨します。これにより、関数型プログラミングの原則に従い、より保守しやすく、テストしやすいコードベースになります。
