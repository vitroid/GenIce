# genice3 テストスイート

## 概要

genice3 のシステマティックなテストスイートです。組み合わせ爆発を避けるため、階層的なテスト戦略を採用しています。

## テスト構造

```
tests/genice3/
├── __init__.py
├── conftest.py          # pytest fixtures
├── test_validation.py   # 検証関数の定義
├── test_core.py         # コア機能のテスト
├── test_integration.py  # 統合テスト
└── README.md           # このファイル
```

## テストの階層

### レベル 1: 基本的な動作確認（必須）

- 例外が発生しない
- 出力が生成される

### レベル 2: 出力形式の検証（必須）

- フォーマットが正しい（gromacs, cif など）
- 必須フィールドが存在する

### レベル 3: 物理的・構造的整合性（重要）

- **アイスルール**: 各水分子から 2 本の水素結合が出て、2 本が入る
- **原子数**: 期待通りの原子数が生成される
- **座標**: 座標がセル内に収まっている

### レベル 4: 再現性（重要）

- 同じ seed で同じ結果が得られる

### レベル 5: 構造的整合性（推奨）

- グラフ構造が正しい
- グラフが連結している

## 実行方法

### 前提条件

**プロジェクトのルートディレクトリ（`GenIce2`）で実行する必要があります。**

```
/Users/matto/Dropbox/gitbox/GenIce2/  ← ここで実行
├── genice3/                          ← genice3モジュール
│   ├── __init__.py
│   ├── genice.py
│   └── ...
└── tests/
    └── genice3/                      ← テストファイル
        ├── conftest.py
        └── test_*.py
```

### 仮想環境のアクティベート

```bash
# プロジェクトのルートディレクトリで
source .venv/bin/activate
# または
. .venv/bin/activate
```

### テストの実行

**推奨方法: プロジェクトのルートディレクトリで実行**

```bash
# プロジェクトのルートディレクトリ（GenIce2）で実行
cd /Users/matto/Dropbox/gitbox/GenIce2

# すべてのテストを実行
pytest tests/genice3/ -v

# 特定のテストファイルを実行
pytest tests/genice3/test_core.py -v

# 特定のテスト関数を実行
pytest tests/genice3/test_core.py::test_genice3_initialization -v
```

**注意**: `pytest.ini`で`pythonpath = .`が設定されているので、通常は追加の設定は不要です。

## テストマトリックス

代表的な組み合わせのみをテストします：

- **unitcell**: `1h`, `A15`
- **molecule**: `4site`, `tip4p`
- **exporter**: `gromacs`, `cif`

これにより、全組み合わせ（数十通り）ではなく、主要な組み合わせ（約 10 通り）のみをテストします。

## 検証関数

`test_validation.py`に以下の検証関数が定義されています：

- `validate_basic_execution()`: レベル 1
- `validate_output_format()`: レベル 2
- `validate_ice_rules()`: レベル 3
- `validate_atom_count()`: レベル 3
- `validate_coordinates_in_cell()`: レベル 3
- `validate_graph_structure()`: レベル 5
- `validate_reproducibility()`: レベル 4
- `validate_comprehensive()`: 総合検証

## トラブルシューティング

### モジュールが見つからないエラー

`ModuleNotFoundError: No module named 'genice3.genice'`が発生する場合：

1. **プロジェクトのルートディレクトリで実行しているか確認**

   ```bash
   pwd
   # /Users/matto/Dropbox/gitbox/GenIce2 である必要がある
   ```

2. **`pytest.ini`が存在するか確認**

   ```bash
   ls pytest.ini
   ```

3. **`genice3`ディレクトリが存在するか確認**

   ```bash
   ls genice3/__init__.py
   ```

4. **pytest のキャッシュをクリア**
   ```bash
   pytest --cache-clear
   ```

### 依存パッケージが見つからないエラー

`ModuleNotFoundError: No module named 'numpy'`などが発生する場合：

1. 仮想環境がアクティブになっているか確認
2. 必要なパッケージがインストールされているか確認：
   ```bash
   pip install -r requirements.txt
   # または
   poetry install
   ```

## 今後の拡張

- 回帰テスト（リファレンス出力との比較）
- パフォーマンステスト
- エッジケースのテスト
- プラグイン固有のテスト
