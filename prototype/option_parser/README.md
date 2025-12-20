# プールベースオプションパーサーのテスト

このディレクトリには、新しいプールベースのオプションパーサーの実装とテストが含まれています。

## 仕様

このパーサーは、OPTION_HANDLING_PLANS.md で提案されている「さらにユーザーの学習コストが低い方法」を実装しています。

### 主な特徴

1. **基底レベルのオプションを処理**: `--rep`, `--seed`, `--depol_loop`, `--spot_anion`, `--spot_cation`など、genice3 が直接認識するオプションを処理します。

2. **未処理オプションの順次処理**: 基底レベルで処理されなかったオプションは、指定されたプラグインに順次渡されます。例えば、unitcell プラグインが指定されている場合は、まず unitcell プラグインに渡され、unitcell プラグインが処理しなかったオプションは、次に exporter プラグインが指定されている場合は exporter プラグインに渡されます。プラグインは指定される場合もされない場合もあり、指定されたプラグインのみが処理に参加します。

3. **YAML サポート**: YAML 設定ファイルからもオプションを読み込めます。

4. **複数回指定のサポート**: 同じオプションが複数回指定された場合（例：`--guest A12=me --guest A14=et`）、リストとして扱います。

## 使い方

### コマンドライン引数からのパース

```python
from pool_parser import PoolBasedParser

parser = PoolBasedParser()
args = [
    "A15",
    "--rep", "2", "2", "2",
    "--exporter", "gromacs",
    "--seed", "42",
    "--depol_loop", "2000",
    "--spot_anion", "1=Cl",
    "--spot_cation", "5=Na",
    "--shift", "0.1", "0.1", "0.1",
    "--anion", "15=Cl",
    "--cation", "21=Na",
    "--density", "0.8",
    "--guest", "A12=me",
    "--guest", "A14=et",
    "--spot_guest", "0=foursite",
    "--water_model", "foursite",
    "--type", "ice",
]
parser.parse_args(args)

result = parser.get_result()
print(result)
```

### YAML ファイルからのパース

```python
from pool_parser import PoolBasedParser

parser = PoolBasedParser()
parser.parse_yaml("test_config.yaml")

result = parser.get_result()
print(result)
```

## テストの実行

```bash
# pytestを使用する場合
pytest test_pool_parser.py -v

# または、直接実行（yamlモジュールが不要なテストのみ）
python3 test_runner.py
```

## パース結果の構造

```python
{
    "base_options": {
        "seed": "42",
        "depol_loop": "2000",
        "replication_factors": ("2", "2", "2"),
        "spot_anion": {"1": "Cl"},
        "spot_cation": {"5": "Na"},
    },
    "unitcell": {
        "name": "A15",
        "options": {
            "shift": ("0.1", "0.1", "0.1"),
            "anion": "15=Cl",
            "cation": "21=Na",
            "density": "0.8",
            # unitcellプラグインが処理しなかったオプションも含まれる
            "guest": ["A12=me", "A14=et"],
            "spot_guest": "0=foursite",
            "water_model": "foursite",
            "type": "ice",
        }
    },
    "exporter": {
        "name": "gromacs",
        "options": {}  # 実際の実装では、unitcellが処理しなかったオプションが入る
    },
    "unprocessed_options": {}  # 最終的に処理されなかったオプション
}
```

## プラグインの実装

実際のプラグインは、`parse_options`関数を実装し、オプションの辞書を受け取って、処理したオプションと処理しなかったオプションを返します。

### A15 プラグインの例

```python
from pool_parser import OptionPool

def parse_options(options: Dict[str, Any]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    pool = OptionPool(options)
    processed = {}

    # shift, density, anion, cation を処理
    if "shift" in pool:
        processed["shift"] = pool.get("shift")
    # ... 他のオプションも同様

    # 処理されなかったオプションを返す
    unprocessed = pool.get_unprocessed()
    return processed, unprocessed
```

### gromacs プラグインの例

`guest`, `spot_guest`, `water_model` オプションを処理します。`type` オプションは `water_model` が `foursite` の場合、`foursite` molecule プラグインで処理されます。

### foursite molecule プラグインの例

`type` オプションを処理します。`gromacs` プラグインから呼び出され、処理したオプションは `water_model` の `options` に統合されます。

### プラグインを使用したテスト

`test_with_plugins.py` を実行すると、実際のプラグインを使用してオプションが処理されることを確認できます:

```bash
python3 test_with_plugins.py
```

このテストでは、unitcell プラグイン（A15）が処理しなかったオプションが exporter プラグイン（gromacs）に渡されることを確認します。

## 注意事項

- この実装は**プロトタイプ**です。実際の genice3 に組み込むには、さらなる改善が必要です。
- プラグインは `OptionPool` を使用してオプションを処理し、処理しなかったオプションを返します。
- `[plugin --option]`形式のパースは部分的に実装されていますが、完全ではありません。
