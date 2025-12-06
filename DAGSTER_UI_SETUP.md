# Dagster UI の起動方法

## 前提条件

1. Dagster の依存関係をインストール:

```bash
poetry install
# または
pip install dagster dagster-webserver
```

## 起動方法

### 方法 1: Dagster CLI を使用（推奨）

```bash
dagster dev -m genice2.genice3_dagster
```

### 方法 2: Python スクリプトを使用

```bash
python run_dagster.py
```

### 方法 3: 直接 Python から起動

```bash
python -m dagster dev -m genice2.genice3_dagster
```

## UI へのアクセス

起動後、ブラウザで以下の URL にアクセス:

```
http://localhost:3000
```

## 使用方法

1. **Asset グラフの表示**: 左側のメニューから「Assets」を選択すると、すべての asset とその依存関係がグラフとして表示されます。

2. **Asset の実行**: 各 asset をクリックして、実行ボタンを押すと、その asset とその依存関係が実行されます。

3. **設定の変更**: 現在の実装では、UI で表示される asset はデフォルト値を使用します。パラメータを変更するには、CLI を使用してください。

## 注意事項

### UI で表示される asset と CLI で使用される asset の違い

- **UI で表示される asset**: `genice3_config`、`unitcell_asset`、`cell`など、静的に定義された asset。これらはデフォルト値を使用します。
- **CLI で使用される asset**: `dynamic_*`で始まる asset。これらは CLI 実行時に動的に作成され、UI には表示されません。CLI パラメータに基づいて設定されます。

UI でパラメータを変更可能にするには、Dagster の Config や Resource を使用する必要があります。現在の実装では、UI ではデフォルト値で実行でき、パラメータを変更するには CLI を使用してください。

## 設定のカスタマイズ

`genice3_config`と`unitcell_asset`の asset は、以下のパラメータを受け取ります:

### genice3_config

- `depol_loop`: 分極ループ回数（デフォルト: 1000）
- `replication_matrix`: 複製行列（3x3 のリスト、デフォルト: 単位行列）

### unitcell_asset

- `unitcell_type`: 単位胞タイプ（"A15"または"Ice1h"、デフォルト: "A15"）
- `shift`: 単位胞のシフト（デフォルト: (0.0, 0.0, 0.0)）
- `assess_cages`: ケージを評価するかどうか（デフォルト: False）
- `anions`: アニオン位置の辞書（デフォルト: {}）
- `cations`: カチオン位置の辞書（デフォルト: {}）
- `density`: 密度（オプション）

## トラブルシューティング

### インポートエラーが発生する場合

Dagster がインストールされていない可能性があります:

```bash
poetry install
```

### ポートが既に使用されている場合

別のポートを指定:

```bash
dagster dev -m genice2.genice3_dagster --port 3001
```

### Asset が実行されない場合

依存関係が正しく定義されているか確認してください。Dagster UI の「Assets」タブで、依存関係グラフを確認できます。
