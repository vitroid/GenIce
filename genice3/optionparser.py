"""
プラグインオプションの共通パーサー

基本的なオプションの解釈（階層的キー名、丸括弧配列、カンマ区切りなど）を提供。
複雑な記述（混合ガスなど）はプラグイン側で解釈する。

使用例:
    >>> from genice3.optionparser import parse_options

    >>> # 基本的なオプションのパース（階層的キー名は自動的に階層的な辞書構造に変換）
    >>> options = parse_options("guest.A12=me,shift=(0.1,0.1,0.1),verbose")
    >>> options
    {'guest': {'A12': 'me'}, 'shift': [0.1, 0.1, 0.1], 'verbose': True}

    >>> # 混合ガスのような複雑な記述は、値として文字列のまま残る
    >>> options = parse_options("guest.A12=me*0.3[monatomic]+et*0.6[molecular]")
    >>> options['guest']['A12']
    'me*0.3[monatomic]+et*0.6[molecular]'  # プラグイン側で解釈
"""

from typing import Dict, List, Any, Union
from logging import getLogger

# from genice3.genice import GuestSpec
# from genice3.plugin import safe_import
# from genice3.molecule import Molecule

logger = getLogger(__name__)


def parse_options(option_string: str) -> Dict[str, Any]:
    """
    プラグインオプション文字列をパースして辞書に変換

    サポートする形式:
    - 階層的キー名: `guest.A12=me` → `{"guest": {"A12": "me"}}`
    - 丸括弧配列: `shift=(0.1,0.1,0.1), rep=(3,3,3)`
    - 単一値: `water_model=tip4p`
    - フラグ: `verbose`、`?`、`help?`、`cage?` など（値なしでTrueが設定される）

    Args:
        option_string: オプション文字列（例: "guest.A12=me,shift=(0.1,0.1,0.1),verbose"）

    Returns:
        パース済みのオプション辞書（階層的キー名は自動的に階層的な辞書構造に変換される）

    Examples:
        >>> parse_options("guest.A12=me,shift=(0.1,0.1,0.1)")
        {"guest": {"A12": "me"}, "shift": [0.1, 0.1, 0.1]}

        >>> parse_options("rep=(3,3,3),verbose")
        {"rep": [3, 3, 3], "verbose": True}

        >>> parse_options("guest.A12=me,guest.A14=et")
        {"guest": {"A12": "me", "A14": "et"}}
    """
    result = {}

    if not option_string or not option_string.strip():
        return result

    option_string = option_string.strip()

    # カンマでオプションを分割（ただし、丸括弧内のカンマは無視）
    options = _split_options(option_string)

    for opt in options:
        opt = opt.strip()
        if not opt:
            continue

        if "=" in opt:
            key, value = opt.split("=", 1)
            key = key.strip()
            value = value.strip()

            # 値をパース
            parsed_value = None
            if value.startswith("(") and value.endswith(")"):
                # 配列/ベクトル: shift=(0.1,0.1,0.1)
                parsed_value = _parse_array(value)
            else:
                # 単一値
                parsed_value = value

            # 階層的キー名を階層的な辞書構造に変換
            _set_nested_value(result, key, parsed_value)
        else:
            # フラグ（値なし）
            _set_nested_value(result, opt, True)

    return result


def _set_nested_value(result: Dict[str, Any], key: str, value: Any) -> None:
    """
    階層的キー名で値を設定（階層的な辞書構造に変換）

    Args:
        result: 結果辞書（更新される）
        key: 階層的キー名（例: "guest.A12" または "spot_guest.0"）
        value: 設定する値

    Examples:
        >>> result = {}
        >>> _set_nested_value(result, "guest.A12", "me")
        >>> result
        {"guest": {"A12": "me"}}

        >>> _set_nested_value(result, "guest.A14", "et")
        >>> result
        {"guest": {"A12": "me", "A14": "et"}}
    """
    if "." in key:
        # 階層的キー名: guest.A12 → {"guest": {"A12": value}}
        parts = key.split(".")
        current = result
        for part in parts[:-1]:
            if part not in current:
                current[part] = {}
            elif not isinstance(current[part], dict):
                # 既存の値が辞書でない場合は上書き
                current[part] = {}
            current = current[part]
        # 最後の部分に値を設定
        current[parts[-1]] = value
    else:
        # 単純なキー名
        result[key] = value


def _split_options(option_string: str) -> List[str]:
    """
    オプション文字列を分割（丸括弧内のカンマは無視）

    例: "guest.A12=me,shift=(0.1,0.1,0.1),rep=(3,3,3)"
    -> ["guest.A12=me", "shift=(0.1,0.1,0.1)", "rep=(3,3,3)"]

    ただし、角括弧内のカンマも無視する必要がある（ネストしたプラグインのオプション）
    例: "guest.A12=me[monatomic,verbose],shift=(0.1,0.1,0.1)"
    -> ["guest.A12=me[monatomic,verbose]", "shift=(0.1,0.1,0.1)"]
    """
    result = []
    current = ""
    depth_paren = 0  # 丸括弧の深さ
    depth_bracket = 0  # 角括弧の深さ

    for char in option_string:
        if char == "(":
            depth_paren += 1
            current += char
        elif char == ")":
            depth_paren -= 1
            current += char
        elif char == "[":
            depth_bracket += 1
            current += char
        elif char == "]":
            depth_bracket -= 1
            current += char
        elif char == "," and depth_paren == 0 and depth_bracket == 0:
            # トップレベルのカンマで分割
            if current:
                result.append(current.strip())
                current = ""
        else:
            current += char

    if current:
        result.append(current.strip())

    return result


def _parse_array(array_string: str) -> List[Union[float, int, str]]:
    """
    丸括弧配列文字列をパース

    Args:
        array_string: 配列文字列（例: "(0.1,0.1,0.1)" または "(3,3,3)"）

    Returns:
        パース済みの配列（数値は自動的に変換）

    Examples:
        >>> _parse_array("(0.1,0.1,0.1)")
        [0.1, 0.1, 0.1]

        >>> _parse_array("(3,3,3)")
        [3, 3, 3]

        >>> _parse_array("(a,b,c)")
        ['a', 'b', 'c']
    """
    if not (array_string.startswith("(") and array_string.endswith(")")):
        raise ValueError(f"Invalid array format: {array_string}")

    content = array_string[1:-1]  # 括弧を除去
    if not content.strip():
        return []

    elements = [elem.strip() for elem in content.split(",")]
    result = []

    for elem in elements:
        if not elem:
            continue
        # 数値に変換できる場合は変換
        try:
            # 浮動小数点数として試す
            result.append(float(elem))
        except ValueError:
            try:
                # 整数として試す
                result.append(int(elem))
            except ValueError:
                # 変換できない場合は文字列のまま
                result.append(elem)

    return result


def _set_nested_value(result: Dict[str, Any], key: str, value: Any) -> None:
    """
    階層的キー名で値を設定（階層的な辞書構造に変換）

    `guest.A12=me` → `{"guest": {"A12": "me"}}`

    Args:
        result: 結果辞書（更新される）
        key: 階層的キー名（例: "guest.A12" または "spot_guest.0"）
        value: 設定する値

    Examples:
        >>> result = {}
        >>> _set_nested_value(result, "guest.A12", "me")
        >>> result
        {"guest": {"A12": "me"}}

        >>> _set_nested_value(result, "guest.A14", "et")
        >>> result
        {"guest": {"A12": "me", "A14": "et"}}

        >>> result = {}
        >>> _set_nested_value(result, "shift", [0.1, 0.1, 0.1])
        >>> result
        {"shift": [0.1, 0.1, 0.1]}
    """
    if "." in key:
        # 階層的キー名: guest.A12 → {"guest": {"A12": value}}
        parts = key.split(".")
        current = result
        for part in parts[:-1]:
            if part not in current:
                current[part] = {}
            elif not isinstance(current[part], dict):
                # 既存の値が辞書でない場合は上書き
                current[part] = {}
            current = current[part]
        # 最後の部分に値を設定
        current[parts[-1]] = value
    else:
        # 単純なキー名
        result[key] = value


def parse_hierarchical_key(key: str) -> List[str]:
    """
    階層的キー名をドットで分割

    Args:
        key: 階層的キー名（例: "guest.A12" または "spot_guest.0"）

    Returns:
        キーの階層リスト

    Examples:
        >>> parse_hierarchical_key("guest.A12")
        ['guest', 'A12']

        >>> parse_hierarchical_key("spot_guest.0")
        ['spot_guest', '0']
    """
    return key.split(".")


def get_nested_value(
    options: Dict[str, Any], key_path: List[str], default: Any = None
) -> Any:
    """
    階層的キーで値を取得

    Args:
        options: オプション辞書
        key_path: キーの階層リスト（例: ["guest", "A12"]）
        default: デフォルト値

    Returns:
        取得した値、またはデフォルト値

    Examples:
        >>> options = {"guest.A12": "me", "guest.A14": "et"}
        >>> get_nested_value(options, ["guest", "A12"])
        'me'

        >>> get_nested_value(options, ["guest", "A15"], "none")
        'none'
    """
    key = ".".join(key_path)
    return options.get(key, default)


def extract_options_with_prefix(options: Dict[str, Any], prefix: str) -> Dict[str, Any]:
    """
    指定されたプレフィックスで始まるオプションを抽出

    Args:
        options: オプション辞書
        prefix: プレフィックス（例: "guest." または "spot_guest."）

    Returns:
        プレフィックスを除去したキーを持つ新しい辞書

    Examples:
        >>> options = {"guest.A12": "me", "guest.A14": "et", "shift": "(0.1,0.1,0.1)"}
        >>> extract_options_with_prefix(options, "guest.")
        {"A12": "me", "A14": "et"}
    """
    result = {}
    for key, value in options.items():
        if key.startswith(prefix):
            new_key = key[len(prefix) :]
            result[new_key] = value
    return result


def group_by_hierarchical_prefix(
    options: Dict[str, Any], prefix: str
) -> Dict[str, Dict[str, Any]]:
    """
    階層的プレフィックスでオプションをグループ化

    Args:
        options: オプション辞書
        prefix: プレフィックス（例: "guest."）

    Returns:
        グループ化された辞書（例: {"A12": {...}, "A14": {...}}）

    Examples:
        >>> options = {"guest.A12": "me", "guest.A12.monatomic": True, "guest.A14": "et"}
        >>> group_by_hierarchical_prefix(options, "guest.")
        {"A12": {"": "me", "monatomic": True}, "A14": {"": "et"}}
    """
    result = {}
    prefix_len = len(prefix)

    for key, value in options.items():
        if not key.startswith(prefix):
            continue

        remaining = key[prefix_len:]
        if "." in remaining:
            # さらに階層がある場合: guest.A12.monatomic
            parts = remaining.split(".", 1)
            group_key = parts[0]
            sub_key = parts[1]

            if group_key not in result:
                result[group_key] = {}
            result[group_key][sub_key] = value
        else:
            # トップレベルの値: guest.A12
            if remaining not in result:
                result[remaining] = {}
            result[remaining][""] = value

    return result
