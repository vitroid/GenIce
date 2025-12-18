# dictparser モジュール

文字列の辞書をパースして、キーに応じて値を整数型やリスト型などに解釈するユーティリティです。

## 概要

genice3 の中で、文字列の辞書（例: `{"width": "100", "height": "200", "HB": "0.4"}`）をパースして、キーに応じて値を型変換する処理が何度も繰り返されています。このモジュールは、そのような処理を簡潔に記述するためのツールです。

## 基本的な使い方

### 方法 1: シンプルな型マッピング（推奨）

最もシンプルで使いやすい方法です。

```python
from genice3.dictparser import parse_dict_options

# 型マッピングを定義
type_map = {
    "width": int,
    "height": int,
    "margin": int,
    "HB": float,
    "O": float,
}

# デフォルト値を定義
defaults = {
    "width": 0,
    "height": 0,
    "margin": 0,
    "HB": 0.4,
    "O": 0.06,
}

# リストとして扱うキーを定義（オプション）
list_keys = {
    "rotate": ",",  # カンマで分割
}

# 文字列の辞書をパース
data = {
    "width": "100",
    "height": "200",
    "margin": "10",
    "HB": "0.5",
    "rotate": "x30,y45,z60",
}

options = parse_dict_options(
    data,
    type_map=type_map,
    defaults=defaults,
    list_keys=list_keys,
)

# 属性としてアクセス可能
print(options.width)    # 100 (int)
print(options.height)    # 200 (int)
print(options.HB)        # 0.5 (float)
print(options.rotate)    # ['x30', 'y45', 'z60'] (list)
print(options.O)         # 0.06 (default値)
```

### 方法 2: デコレータを使う方法

より宣言的なアプローチです。

```python
from genice3.dictparser import DictParser

parser = DictParser()

@parser.option("width", type=int, default=0)
@parser.option("height", type=int, default=0)
@parser.option("margin", type=int, default=0)
@parser.option("HB", type=float, default=0.4)
@parser.option("O", type=float, default=0.06)
@parser.option("rotate", multiple=True, separator=",")
def parse_options(**kwargs):
    return parser.parse(kwargs)

# 使用
data = {
    "width": "100",
    "height": "200",
    "margin": "10",
    "HB": "0.5",
    "rotate": "x30,y45,z60",
}

options = parse_options(**data)

print(options.width)    # 100 (int)
print(options.height)    # 200 (int)
print(options.HB)        # 0.5 (float)
print(options.rotate)    # ['x30', 'y45', 'z60'] (list)
```

## 実際の使用例

### svg.py の parse_options 関数をリファクタリングする場合

```python
from genice3.dictparser import parse_dict_options

def parse_options(**kwargs):
    logger = getLogger()

    # 型マッピングを定義
    type_map = {
        "width": int,
        "height": int,
        "margin": int,
        "HB": float,
        "O": float,
        "OH": float,
        "H": float,
        "encode": lambda x: bool(x) if isinstance(x, str) else bool(x),
    }

    # デフォルト値を定義
    defaults = {
        "encode": True,
        "poly": False,
        "shadow": None,
        "oxygen": 0.06,
        "HB": 0.4,
        "OH": 0.5,
        "hydrogen": 0,
        "arrows": False,
        "bgcolor": None,
        "width": 0,
        "height": 0,
        "margin": 0,
    }

    # 基本的な型変換を適用
    options = parse_dict_options(kwargs, type_map=type_map, defaults=defaults)

    # 複雑な処理（回転行列など）は個別に処理
    proj = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])

    for key, value in kwargs.items():
        if key == "rotate":
            # 複雑な処理は従来通り
            values = value.split(",")
            for value in values:
                # ... 回転行列の計算 ...
                pass
        elif key == "polygon" and value is True:
            options.options["poly"] = True
        elif key == "arrows" and value is True:
            options.options["arrows"] = True

    # Optionsオブジェクトを作成
    return Options(
        encode=options.get("encode", True),
        poly=options.get("poly", False),
        shadow=options.get("shadow"),
        oxygen=options.get("O", 0.06),
        HB=options.get("HB", 0.4),
        OH=options.get("OH", 0.5),
        hydrogen=options.get("H", 0),
        arrows=options.get("arrows", False),
        bgcolor=options.get("bg"),
        proj=proj,
        width=options.width - options.margin * 2,
        height=options.height - options.margin * 2,
        margin=options.margin,
        unprocessed=options.unprocessed,
    )
```

## click を利用する方法について

click は主にコマンドライン引数のパースに使われます。辞書のパースにも応用できますが、以下の理由から推奨されません：

1. click はコマンドライン引数用に設計されている
2. 辞書を直接パースするには、click.Context を使う必要があり、複雑になる
3. dictparser モジュールの方が、辞書のパースに特化しており、シンプルで使いやすい

ただし、コマンドライン引数と辞書のパースを統一したい場合は、click の型変換関数を再利用できます：

```python
import click
from genice3.dictparser import parse_dict_options

# clickの型変換関数を再利用
type_map = {
    "width": click.INT,
    "height": click.INT,
    "HB": click.FLOAT,
}

options = parse_dict_options(data, type_map=type_map, defaults=defaults)
```

## API リファレンス

### `parse_dict_options(data, type_map=None, defaults=None, list_keys=None)`

シンプルな方法で辞書をパースします。

**引数:**

- `data`: 文字列の辞書
- `type_map`: キー名から型変換関数へのマッピング（例: `{"width": int, "HB": float}`）
- `defaults`: デフォルト値の辞書
- `list_keys`: リストとして扱うキーとその区切り文字のマッピング（例: `{"rotate": ","}`）

**戻り値:**

- `ParsedOptions`オブジェクト（属性としてアクセス可能）

### `DictParser`

デコレータベースのパーサー。

**メソッド:**

- `option(name, type=str, default=None, required=False, help=None, multiple=False, separator=",")`: オプションを登録するデコレータ
- `parse(data)`: 辞書をパースする

### `ParsedOptions`

パース済みのオプションオブジェクト。

**属性:**

- `options`: パース済みのオプション辞書
- `unprocessed`: 未処理のオプション辞書

**メソッド:**

- `get(name, default=None)`: 辞書のように値を取得
- 属性として直接アクセス可能（例: `options.width`）
