"""
文字列の辞書をパースして、キーに応じて値を型変換するデコレータ

clickに似たデコレータベースのアプローチで、文字列の辞書をパースします。

基本的な使用例:
    >>> from genice3.dictparser import parse_dict_options

    >>> # シンプルな方法：型マッピング辞書を定義
    >>> type_map = {
    ...     "width": int,
    ...     "height": int,
    ...     "HB": float,
    ...     "rotate": lambda x: x.split(","),
    ... }
    >>> defaults = {"width": 0, "height": 0, "HB": 0.4}
    >>>
    >>> options = parse_dict_options(
    ...     {"width": "100", "height": "200", "HB": "0.5", "rotate": "x30,y45"},
    ...     type_map=type_map,
    ...     defaults=defaults
    ... )
    >>> options.width
    100
    >>> options.height
    200
    >>> options.HB
    0.5
    >>> options.rotate
    ['x30', 'y45']

デコレータを使った方法:
    >>> from genice3.dictparser import DictParser

    >>> parser = DictParser()
    >>>
    >>> @parser.option("width", type=int, default=0)
    ... @parser.option("height", type=int, default=0)
    ... @parser.option("HB", type=float, default=0.4)
    ... @parser.option("rotate", multiple=True, separator=",")
    ... def parse_options(**kwargs):
    ...     return parser.parse(kwargs)
    ...
    >>> options = parse_options(width="100", height="200", HB="0.5", rotate="x30,y45")
    >>> options.width
    100
    >>> options.rotate
    ['x30', 'y45']
"""

from typing import Any, Callable, Dict, Optional, Type, Union
from dataclasses import dataclass, field
from functools import wraps
from logging import getLogger

logger = getLogger(__name__)


@dataclass
class Option:
    """オプションの定義"""

    name: str
    type: Callable = str
    default: Any = None
    required: bool = False
    help: Optional[str] = None
    multiple: bool = False  # リストとして扱うか
    separator: str = ","  # リストの区切り文字


class DictParser:
    """辞書パーサーのクラス"""

    def __init__(self):
        self.options: Dict[str, Option] = {}

    def option(
        self,
        name: str,
        type: Callable = str,
        default: Any = None,
        required: bool = False,
        help: Optional[str] = None,
        multiple: bool = False,
        separator: str = ",",
    ):
        """
        オプションを登録するデコレータ

        Args:
            name: オプション名（辞書のキー）
            type: 型変換関数（int, float, str, またはカスタム関数）
            default: デフォルト値
            required: 必須かどうか
            help: ヘルプ文字列
            multiple: リストとして扱うか（separatorで分割）
            separator: リストの区切り文字（multiple=Trueの場合）
        """

        def decorator(func):
            self.options[name] = Option(
                name=name,
                type=type,
                default=default,
                required=required,
                help=help,
                multiple=multiple,
                separator=separator,
            )
            return func

        return decorator

    def parse(self, data: Dict[str, Any]) -> Any:
        """
        辞書をパースして、型変換されたオプションオブジェクトを返す

        Args:
            data: 文字列の辞書（例: {"width": "100", "height": "200"}）

        Returns:
            パース済みのオプションオブジェクト（属性としてアクセス可能）
        """
        result = {}
        unprocessed = {}

        # 登録されたオプションを処理
        for name, option in self.options.items():
            if name in data:
                value = data[name]
                try:
                    if option.multiple:
                        # リストとして扱う
                        if isinstance(value, str):
                            # 空白で分割する場合（Noneまたは" "）は、split()（引数なし）を使用
                            if option.separator is None or option.separator == " ":
                                parts = value.split()
                            else:
                                parts = value.split(option.separator)
                            result[name] = [
                                option.type(v.strip()) for v in parts if v.strip()
                            ]
                        elif isinstance(value, (list, tuple)):
                            result[name] = [option.type(v) for v in value]
                        else:
                            result[name] = [option.type(value)]
                    else:
                        # 単一値として扱う
                        if value is True or value is False:
                            # フラグ（値なしでTrueが設定される場合）
                            result[name] = value
                        else:
                            result[name] = option.type(value)
                except (ValueError, TypeError) as e:
                    logger.warning(
                        f"Failed to convert {name}={value} to {option.type}: {e}"
                    )
                    result[name] = (
                        option.default if option.default is not None else value
                    )
            elif option.required:
                raise ValueError(f"Required option '{name}' is missing")
            elif option.default is not None:
                result[name] = option.default

        # 未処理のオプションを保持
        for key, value in data.items():
            if key not in self.options:
                unprocessed[key] = value

        # オプションオブジェクトを作成
        return ParsedOptions(result, unprocessed)


@dataclass
class ParsedOptions:
    """パース済みのオプションオブジェクト"""

    options: Dict[str, Any] = field(default_factory=dict)
    unprocessed: Dict[str, Any] = field(default_factory=dict)

    def __getattr__(self, name: str) -> Any:
        """属性としてアクセス可能にする"""
        if name in self.options:
            return self.options[name]
        elif name in self.unprocessed:
            return self.unprocessed[name]
        else:
            raise AttributeError(f"Option '{name}' not found")

    def __contains__(self, name: str) -> bool:
        """in演算子をサポート"""
        return name in self.options or name in self.unprocessed

    def get(self, name: str, default: Any = None) -> Any:
        """辞書のように値を取得"""
        if name in self.options:
            return self.options[name]
        elif name in self.unprocessed:
            return self.unprocessed[name]
        else:
            return default


def parse_dict_options(
    data: Dict[str, Any],
    type_map: Optional[Dict[str, Callable]] = None,
    defaults: Optional[Dict[str, Any]] = None,
    list_keys: Optional[Dict[str, str]] = None,
) -> ParsedOptions:
    """
    シンプルな方法：型マッピング辞書を使って辞書をパース

    Args:
        data: 文字列の辞書
        type_map: キー名から型変換関数へのマッピング（例: {"width": int, "HB": float}）
        defaults: デフォルト値の辞書
        list_keys: リストとして扱うキーとその区切り文字のマッピング（例: {"rotate": ","}）

    Returns:
        パース済みのオプションオブジェクト

    使用例:
        >>> type_map = {"width": int, "height": int, "HB": float}
        >>> defaults = {"width": 0, "height": 0}
        >>> options = parse_dict_options(
        ...     {"width": "100", "height": "200", "HB": "0.5"},
        ...     type_map=type_map,
        ...     defaults=defaults
        ... )
        >>> options.width
        100
    """
    type_map = type_map or {}
    defaults = defaults or {}
    list_keys = list_keys or {}

    result = {}
    unprocessed = {}

    for key, value in data.items():
        if key in list_keys:
            # リストとして扱う（type_mapに含まれていなくても処理）
            separator = list_keys[key]
            type_func = type_map.get(key, str)  # type_mapにない場合はstrを使用
            try:
                if isinstance(value, str):
                    # 空白で分割する場合（Noneまたは" "）は、split()（引数なし）を使用
                    # これにより、任意の空白文字（スペース、タブ、改行など）で分割し、
                    # 連続する空白も1つの区切りとして扱う
                    if separator is None or separator == " ":
                        parts = value.split()
                    else:
                        parts = value.split(separator)
                    result[key] = [type_func(v.strip()) for v in parts if v.strip()]
                elif isinstance(value, (list, tuple)):
                    result[key] = [type_func(v) for v in value]
                else:
                    result[key] = [type_func(value)]
            except (ValueError, TypeError) as e:
                logger.warning(f"Failed to convert {key}={value} to {type_func}: {e}")
                result[key] = defaults.get(key, value)
        elif key in type_map:
            # 型変換を適用（リストではない）
            try:
                # 単一値として扱う
                if value is True or value is False:
                    result[key] = value
                else:
                    result[key] = type_map[key](value)
            except (ValueError, TypeError) as e:
                logger.warning(
                    f"Failed to convert {key}={value} to {type_map[key]}: {e}"
                )
                result[key] = defaults.get(key, value)
        elif key in defaults:
            # デフォルト値のみ設定
            result[key] = defaults[key]
        else:
            # 未処理のオプション
            unprocessed[key] = value

    # デフォルト値で未設定のキーを追加
    for key, default_value in defaults.items():
        if key not in result:
            result[key] = default_value

    return ParsedOptions(result, unprocessed)
