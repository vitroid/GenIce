# プールベースオプションパーサーのテスト

このディレクトリには、新しいプールベースのオプションパーサーの実装とテストが含まれています。

## 仕様

このパーサーは、OPTION_HANDLING_PLANS.md で提案されている「さらにユーザーの学習コストが低い方法」を実装しています。

### 主な特徴

1. **基底レベルのオプションを処理**: `--rep`, `--seed`, `--depol_loop`, `--spot_anion`, `--spot_cation`など、genice3 が直接認識するオプションを処理します。

2. **未処理オプションの順次処理**: 基底レベルで処理されなかったオプションは、unitcell プラグインに渡されます。unitcell プラグインが処理しなかったオプションは、exporter プラグインに渡されます。

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
    "--spot_guest", "0=4site",
    "--water_model", "4site",
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
            "spot_guest": "0=4site",
            "water_model": "4site",
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

## 注意事項

- この実装は**プロトタイプ**です。実際の genice3 に組み込むには、さらなる改善が必要です。
- 現在の実装では、unitcell プラグインが処理しなかったオプションを自動的に exporter プラグインに渡すロジックは含まれていません。実際のプラグイン処理と統合する必要があります。
- `[plugin --option]`形式のパースは部分的に実装されていますが、完全ではありません。
