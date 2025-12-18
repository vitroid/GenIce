"""
dictparserの使用例

このファイルは、dictparserモジュールの使い方を示す例です。
実際のコードでは、このファイルを削除して構いません。
"""

from genice3.dictparser import DictParser, parse_dict_options


# 方法1: シンプルな型マッピングを使う方法（推奨）
def example_simple():
    """シンプルな方法の例"""
    # 型マッピングを定義
    type_map = {
        "width": int,
        "height": int,
        "margin": int,
        "HB": float,
        "O": float,
        "OH": float,
        "H": float,
    }

    # デフォルト値を定義
    defaults = {
        "width": 0,
        "height": 0,
        "margin": 0,
        "HB": 0.4,
        "O": 0.06,
        "OH": 0.5,
        "H": 0,
    }

    # リストとして扱うキーを定義
    list_keys = {
        "rotate": ",",  # カンマで分割
    }

    # 文字列の辞書をパース
    data = {
        "width": "100",
        "height": "200",
        "margin": "10",
        "HB": "0.5",
        "O": "0.08",
        "rotate": "x30,y45,z60",
    }

    options = parse_dict_options(
        data,
        type_map=type_map,
        defaults=defaults,
        list_keys=list_keys,
    )

    print(f"width: {options.width} (type: {type(options.width)})")
    print(f"height: {options.height} (type: {type(options.height)})")
    print(f"HB: {options.HB} (type: {type(options.HB)})")
    print(f"rotate: {options.rotate} (type: {type(options.rotate)})")
    print(f"O: {options.O} (default from defaults)")
    print(f"unprocessed keys: {options.unprocessed}")


# 方法2: デコレータを使う方法
def example_decorator():
    """デコレータを使う方法の例"""
    parser = DictParser()

    @parser.option("width", type=int, default=0)
    @parser.option("height", type=int, default=0)
    @parser.option("margin", type=int, default=0)
    @parser.option("HB", type=float, default=0.4)
    @parser.option("O", type=float, default=0.06)
    @parser.option("OH", type=float, default=0.5)
    @parser.option("H", type=float, default=0)
    @parser.option("rotate", multiple=True, separator=",")
    def parse_options(**kwargs):
        return parser.parse(kwargs)

    # 使用
    data = {
        "width": "100",
        "height": "200",
        "margin": "10",
        "HB": "0.5",
        "O": "0.08",
        "rotate": "x30,y45,z60",
    }

    options = parse_options(**data)

    print(f"width: {options.width} (type: {type(options.width)})")
    print(f"height: {options.height} (type: {type(options.height)})")
    print(f"HB: {options.HB} (type: {type(options.HB)})")
    print(f"rotate: {options.rotate} (type: {type(options.rotate)})")
    print(f"O: {options.O} (default)")
    print(f"unprocessed keys: {options.unprocessed}")


# 方法3: clickを利用する方法（オプション）
def example_with_click():
    """
    clickを利用する方法の例

    clickは主にコマンドライン引数のパースに使われますが、
    辞書のパースにも応用できます。
    """
    try:
        import click

        # clickの型変換関数を利用
        @click.command()
        @click.option("--width", type=int, default=0)
        @click.option("--height", type=int, default=0)
        @click.option("--HB", type=float, default=0.4)
        def parse_with_click(width, height, HB):
            return {"width": width, "height": height, "HB": HB}

        # ただし、clickはコマンドライン引数用なので、
        # 辞書を直接パースするには、click.Contextを使う必要があります
        # この方法は推奨されません

        print("clickは主にコマンドライン引数のパースに使われます")
        print("辞書のパースには、dictparserモジュールを使用することを推奨します")

    except ImportError:
        print("clickがインストールされていません")


if __name__ == "__main__":
    print("=== 方法1: シンプルな型マッピング ===")
    example_simple()

    print("\n=== 方法2: デコレータ ===")
    example_decorator()

    print("\n=== 方法3: click（参考） ===")
    example_with_click()

