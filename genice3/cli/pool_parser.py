"""
プラグインオプションの共通パーサー

PoolBasedParserを使用した動的プラグインチェーン実行システムを提供。
"""

from typing import Dict, List, Any, Union, Optional, Set, Tuple, Callable
from logging import getLogger, DEBUG, INFO
from pathlib import Path

# from genice3.genice import GuestSpec
# from genice3.plugin import safe_import
# from genice3.molecule import Molecule

logger = getLogger(__name__)

# プラグインのparse_options関数のキャッシュ
_plugin_cache: Dict[Tuple[str, str], Callable] = {}

# YAMLライブラリのインポート
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


# ============================================================================
# OptionPoolクラス定義
# ============================================================================


class OptionPool:
    """未処理のオプションを保持するプール"""

    def __init__(self, options: Dict[str, Any]):
        """
        Args:
            options: オプションの辞書（キーは`--`なしのオプション名、値はその値）
        """
        self.options = options.copy()
        self.processed_keys: Set[str] = set()

    def get(self, key: str, default: Any = None) -> Any:
        """オプションを取得し、処理済みマークを付ける"""
        if key in self.options:
            value = self.options[key]
            self.processed_keys.add(key)
            return value
        return default

    def remove(self, key: str) -> None:
        """オプションを削除"""
        if key in self.options:
            del self.options[key]
            self.processed_keys.add(key)

    def get_unprocessed(self) -> Dict[str, Any]:
        """処理されていないオプションを返す"""
        return {k: v for k, v in self.options.items() if k not in self.processed_keys}

    def has_unprocessed(self) -> bool:
        """処理されていないオプションがあるかどうか"""
        return len(self.get_unprocessed()) > 0

    def __contains__(self, key: str) -> bool:
        return key in self.options


# ============================================================================
# 基底レベルのオプション定義
# ============================================================================

# 基底レベル（genice3）が認識するオプション
BASE_LEVEL_OPTIONS: Set[str] = {
    "rep",
    "replication_factors",
    "replication_matrix",
    "seed",
    "depol_loop",
    "assess_cages",
    "debug",
    "spot_anion",
    "spot_cation",
    "config",
    "exporter",  # exporter名も基底レベルで処理
    "unitcell",  # unitcell名も基底レベルで処理
}


# ============================================================================
# コマンドライン引数と設定ファイルのパース関数（PoolBasedParserに依存しない）
# ============================================================================


def load_config_file(yaml_path: str) -> Dict[str, Any]:
    """
    設定ファイルを読み込んでそのまま返す

    Args:
        yaml_path: YAMLファイルのパス

    Returns:
        YAMLファイルの内容をそのまま返す（YAMLの構造がそのまま内部データ構造になる）
    """
    if not YAML_AVAILABLE:
        return {}

    config_file = Path(yaml_path)
    if not config_file.exists():
        return {}

    with open(config_file, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f) or {}

    return config


def parse_option_string(option_str: str) -> Dict[str, Any]:
    """
    --option value --option2 value2形式の文字列をパース

    Args:
        option_str: オプション文字列

    Returns:
        オプションの辞書（同じオプションが複数回指定された場合はリストとして扱う）
    """
    options = {}
    parts = option_str.split()
    i = 0
    while i < len(parts):
        part = parts[i]
        if part.startswith("--"):
            key = part[2:]
            i += 1
            # 次の値（--で始まっていないもの）を取得
            values = []
            while i < len(parts) and not parts[i].startswith("--"):
                values.append(parts[i])
                i += 1
            # 値が1つの場合はそのまま、複数の場合はタプル
            value = values[0] if len(values) == 1 else tuple(values)

            # 同じオプションが複数回指定された場合はリストに変換
            # （process_option_argと同じロジック）
            if key in options:
                existing = options[key]
                if not isinstance(existing, list):
                    options[key] = [existing]
                options[key].append(value)
            else:
                options[key] = value
        else:
            i += 1
    return options


def parse_bracketed_plugin(arg: str) -> Tuple[str, Dict[str, Any]]:
    """
    [plugin --option value]形式をパース

    Args:
        arg: [plugin --option value]形式の文字列

    Returns:
        (plugin_name, options_dict)
    """
    # [と]を削除
    content = arg[1:-1].strip()
    parts = content.split(None, 1)  # 最初の空白で分割
    plugin_name = parts[0]

    if len(parts) > 1:
        options_str = parts[1]
        options = parse_option_string(options_str)
    else:
        options = {}

    return plugin_name, options


# 短縮形オプションから長い形式へのマッピング
SHORT_OPTION_MAP: Dict[str, str] = {
    "-D": "debug",
    "-s": "seed",
    "-A": "assess_cages",
    "-a": "spot_anion",
    "-c": "spot_cation",
    "-C": "config",
    "-e": "exporter",
    "-h": "help",
    "-V": "version",
}

# フラグ型のオプション（値を持たない）
FLAG_OPTIONS: Set[str] = {"debug"}


def parse_command_line_to_dict(args: List[str]) -> Tuple[Dict[str, Any], Set[str]]:
    """
    コマンドライン引数をパースしてフラットな辞書に変換

    Args:
        args: コマンドライン引数のリスト（--configは除外済みであること）

    Returns:
        (cmdline_dict, cmdline_specified_keys)のタプル
        - cmdline_dict: パース済みのフラットなオプション辞書
        - cmdline_specified_keys: コマンドラインで指定されたキーのセット
    """
    cmdline_dict: Dict[str, Any] = {}
    cmdline_specified_keys: Set[str] = set()
    unitcell_arg = None
    i = 0

    while i < len(args):
        arg = args[i]

        # --configオプションは既に処理済みなのでスキップ
        if arg == "--config" or arg == "-C":
            i += 2  # --config とその値をスキップ
            continue

        # unitcell名の位置引数を処理（最初の位置引数のみ）
        # bracket形式などはパースせず、そのまま辞書に入れる
        # - または -- で始まる場合は位置引数ではない
        if unitcell_arg is None and not arg.startswith("-"):
            cmdline_dict["unitcell"] = arg
            unitcell_arg = arg
            i += 1
            continue

        # --で始まるオプション（長い形式）を処理
        if arg.startswith("--"):
            key = arg[2:]  # --を削除
            cmdline_specified_keys.add(key)
            i += 1

            # 値を取得（次の引数が-で始まっていなければ値として扱う）
            values: List[Any] = []
            while i < len(args) and not args[i].startswith("-"):
                values.append(args[i])
                i += 1

            # 値が1つの場合はそのまま、複数の場合はタプル
            # フラグ型オプションの場合はTrue、それ以外で値がない場合はNone
            if len(values) == 0:
                if key in FLAG_OPTIONS:
                    value = True  # フラグ型オプション
                else:
                    value = None  # 値が必要だが指定されていない
            elif len(values) == 1:
                value = values[0]
            else:
                value = tuple(values)

            # 同じオプションが複数回指定された場合はリストに変換
            if key in cmdline_dict:
                existing = cmdline_dict[key]
                if not isinstance(existing, list):
                    cmdline_dict[key] = [existing]
                cmdline_dict[key].append(value)
            else:
                cmdline_dict[key] = value
        # -で始まるオプション（短縮形式）を処理
        elif arg.startswith("-") and arg != "-":
            # 短縮形オプションを長い形式に変換
            short_arg = arg
            if short_arg in SHORT_OPTION_MAP:
                key = SHORT_OPTION_MAP[short_arg]
                cmdline_specified_keys.add(key)
                i += 1

                # 値を取得（次の引数が-で始まっていなければ値として扱う）
                # フラグ型のオプション（-D, -Aなど）は値なし
                values: List[Any] = []
                while i < len(args) and not args[i].startswith("-"):
                    values.append(args[i])
                    i += 1

                # 値が1つの場合はそのまま、複数の場合はタプル
                # フラグ型オプションの場合はTrue、それ以外で値がない場合はNone
                if len(values) == 0:
                    if key in FLAG_OPTIONS:
                        value = True  # フラグ型オプション
                    else:
                        value = None  # 値が必要だが指定されていない
                elif len(values) == 1:
                    value = values[0]
                else:
                    value = tuple(values)

                # 同じオプションが複数回指定された場合はリストに変換
                if key in cmdline_dict:
                    existing = cmdline_dict[key]
                    if not isinstance(existing, list):
                        cmdline_dict[key] = [existing]
                    cmdline_dict[key].append(value)
                else:
                    cmdline_dict[key] = value
            else:
                # 認識されない短縮形オプションはスキップ
                i += 1
        else:
            i += 1

    return cmdline_dict, cmdline_specified_keys


# ============================================================================
# PoolBasedParser関連（動的プラグインチェーン実行システム）
# ============================================================================


# オプション値の型ごとの処理関数


def process_flag_option(value: Any) -> bool:
    """
    1. 値なし(フラグ)の処理

    Args:
        value: オプションの値（通常はNoneまたは空文字列）

    Returns:
        真偽値（指定されていればTrue）
    """
    return True


def process_string_option(value: Any) -> Any:
    """
    2. 文字列オプションの処理

    Args:
        value: オプションの値（文字列、数値文字列、ブール値文字列など）

    Returns:
        そのままの値（文字列または辞書）
    """
    return value


def process_tuple_option(value: Any) -> list:
    """
    4. リストまたはタプルオプションの処理

    Args:
        value: オプションの値（タプル、リスト、または単一の文字列）

    Returns:
        リストに変換した値
    """
    if isinstance(value, (tuple, list)):
        return list(value)
    else:
        return [value]


def process_keyvalue_option(value: Any, allow_dict: bool = True) -> Dict[str, str]:
    """
    5. "a=b"タイプのオプションの処理

    Args:
        value: オプションの値（文字列、リスト、または辞書）
        allow_dict: Trueの場合、辞書形式も受け入れる（YAMLから来た場合）

    Returns:
        辞書（キーと値のペア）
    """
    result = {}

    if isinstance(value, list):
        # 複数回指定された場合（リスト）- コマンドラインで指定された値
        for item in value:
            if isinstance(item, str) and "=" in item:
                k, v = item.split("=", 1)
                result[k] = v
    elif isinstance(value, str) and "=" in value:
        # 単一のkey=value形式（コマンドライン）
        k, v = value.split("=", 1)
        result[k] = v
    elif isinstance(value, dict) and allow_dict:
        # 辞書形式（設定ファイルから来た値）
        result = value.copy()

    return result


# オプションの型定義用の定数
OPTION_TYPE_FLAG = "flag"  # 1. 値なし(フラグ)
OPTION_TYPE_STRING = "string"  # 2. 文字列
OPTION_TYPE_TUPLE = "tuple"  # 4. リストまたはタプル
OPTION_TYPE_KEYVALUE = "keyvalue"  # 5. "a=b"タイプ


def parse_options_generic(
    options: Dict[str, Any],
    option_specs: Dict[str, str],
    post_processors: Optional[Dict[str, Callable[[Any], Any]]] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    一般的なオプションパーサー

    Args:
        options: オプションの辞書（設定ファイルの値が初期値として含まれる）
        option_specs: オプションの型定義（例：{"shift": OPTION_TYPE_TUPLE, "density": OPTION_TYPE_STRING}）
        post_processors: オプションごとの後処理関数（例：{"shift": lambda x: [float(v) for v in x]}）

    Returns:
        (処理したオプション, 処理しなかったオプション)
    """
    pool = OptionPool(options)
    processed = {}
    post_processors = post_processors or {}

    # option_specsに定義されたオプションを処理
    for key, option_type in option_specs.items():
        if key not in pool:
            continue

        option_value = pool.get(key)

        if option_type == OPTION_TYPE_FLAG:
            processed[key] = process_flag_option(option_value)
        elif option_type == OPTION_TYPE_STRING:
            processed[key] = process_string_option(option_value)
        elif option_type == OPTION_TYPE_TUPLE:
            processed[key] = process_tuple_option(option_value)
        elif option_type == OPTION_TYPE_KEYVALUE:
            processed[key] = process_keyvalue_option(option_value)

        # 後処理がある場合は適用
        if key in post_processors:
            try:
                processed[key] = post_processors[key](processed[key])
            except (ValueError, TypeError):
                pass  # 後処理が失敗しても元の値を保持

    # 処理されなかったオプションを取得
    unprocessed = pool.get_unprocessed()

    return processed, unprocessed


# ============================================================================
# プラグイン読み込み関数
# ============================================================================


def load_plugin(category: str, name: str) -> Optional[Callable]:
    """
    プラグインを動的に読み込む

    Args:
        category: プラグインのカテゴリ（"unitcell", "exporter", "molecule"）
        name: プラグイン名

    Returns:
        parse_options関数、またはNone（読み込めない場合）
    """
    from genice3.plugin import safe_import

    cache_key = (category, name)
    if cache_key in _plugin_cache:
        return _plugin_cache[cache_key]

    try:
        # safe_importを使ってプラグインモジュールを読み込む
        module = safe_import(category, name)
        if hasattr(module, "parse_options"):
            # プラグイン固有のparse_options関数がある場合
            parse_func = module.parse_options
            _plugin_cache[cache_key] = parse_func
            return parse_func
        elif category == "unitcell":
            # unitcellプラグインの場合、基底クラスのparse_optionsを使用
            from genice3.unitcell import UnitCell

            parse_func = UnitCell.parse_options
            _plugin_cache[cache_key] = parse_func
            return parse_func
        elif category == "molecule":
            # moleculeプラグインの場合、基底クラスのparse_optionsを使用
            from genice3.molecule import Molecule

            parse_func = Molecule.parse_options
            _plugin_cache[cache_key] = parse_func
            return parse_func
        elif category == "base":
            # baseオプションの場合、CLI固有のparse_base_optionsを使用
            from genice3.cli.genice import parse_base_options

            parse_func = parse_base_options
            _plugin_cache[cache_key] = parse_func
            return parse_func
    except (ImportError, AttributeError, ModuleNotFoundError, AssertionError):
        pass

    return None


# ============================================================================
# オプション辞書構築関数
# ============================================================================


def build_options_dict(
    config_dict: Dict[str, Any],
    cmdline_dict: Dict[str, Any],
    cmdline_specified_keys: Set[str],
) -> Dict[str, Any]:
    """
    YAML設定とコマンドライン引数から階層的なオプション辞書を構築

    Args:
        config_dict: YAMLファイルから読み込んだ辞書（YAMLの構造をそのまま保持）
        cmdline_dict: コマンドライン引数からパースしたフラットな辞書
        cmdline_specified_keys: コマンドラインで指定されたキーのセット

    Returns:
        階層的なオプション辞書（YAMLの構造に対応）:
        {
            "genice3": {...},  # genice3セクションのオプション
            "unitcell": {...},  # unitcellセクション（nameを含むすべてのオプション）
            "exporter": {...}   # exporterセクション（nameを含むすべてのオプション）
        }
    """
    from genice3.cli.genice import parse_base_options

    # 1. genice3オプションを構築（YAMLから読み込んだgenice3オプションにコマンドラインをマージ）
    genice3_options = config_dict.get("genice3", {}).copy()

    # コマンドラインで指定されたgenice3レベルのオプションで上書き
    for key in BASE_LEVEL_OPTIONS:
        if key in cmdline_dict:
            if key == "rep":
                # repはreplication_factorsのエイリアス
                genice3_options["replication_factors"] = cmdline_dict[key]
            elif key not in ("config", "exporter", "unitcell"):
                genice3_options[key] = cmdline_dict[key]

    # parse_base_optionsでgenice3レベルオプションを処理（統合された辞書を返す）
    genice3_options = parse_base_options(genice3_options)

    # 2. プラグインのオプションを構築
    plugin_types = ["unitcell", "exporter"]
    options_dict: Dict[str, Any] = {
        "genice3": genice3_options,
    }

    # BASE_LEVEL_OPTIONS以外のコマンドライン引数オプションを取得
    # （これらはプラグインのオプションとして扱われる）
    plugin_cmdline_options = {
        k: v
        for k, v in cmdline_dict.items()
        if k not in BASE_LEVEL_OPTIONS and k not in plugin_types
    }

    for plugin_type in plugin_types:
        # YAMLから読み込んだ設定を初期値として使用（構造をそのまま保持）
        plugin_config = config_dict.get(plugin_type, {}).copy()

        # コマンドラインで指定されたプラグイン名で上書き
        if plugin_type in cmdline_dict:
            plugin_spec = cmdline_dict[plugin_type]
            if isinstance(plugin_spec, str):
                if plugin_spec.startswith("[") and plugin_spec.endswith("]"):
                    plugin_name, bracket_options = parse_bracketed_plugin(plugin_spec)
                    # bracket形式: nameを設定し、オプションをマージ
                    plugin_config["name"] = plugin_name
                    plugin_config.update(bracket_options)
                else:
                    # 単純な文字列: nameとして設定
                    plugin_config["name"] = plugin_spec

        # コマンドライン引数のプラグインオプションを追加
        # （どのオプションがどのプラグインに属するかは後でparse_optionsで判定される）
        plugin_config.update(plugin_cmdline_options)

        options_dict[plugin_type] = plugin_config

    logger.debug(f"options_dict: {options_dict}")

    return options_dict


# ============================================================================
# プラグイン実行関数
# ============================================================================


def execute_plugin(
    category: str,
    plugin_name: Optional[str],
    current_options: Dict[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, Any], bool]:
    """
    プラグインを実行して、処理済み/未処理のオプションを返す

    Args:
        category: プラグインのカテゴリ（"unitcell", "exporter", "molecule"など）
        plugin_name: プラグイン名（Noneの場合は実行しない）
        current_options: 現在のオプション辞書

    Returns:
        (processed, unprocessed, executed) のタプル
        - processed: 処理されたオプション
        - unprocessed: 処理されなかったオプション（次のプラグインに渡す）
        - executed: プラグインが実行されたかどうか
    """
    logger = getLogger()
    if not plugin_name:
        logger.debug(f"no plugin name {plugin_name}")
        return {}, current_options, False

    parse_func = load_plugin(category, plugin_name)
    if not parse_func:
        logger.debug(f"Failed to load the plugin {plugin_name}")
        return {}, current_options, False

    processed, unprocessed = parse_func(current_options)

    # debugモードで結果を表示（loggerがDEBUGレベルかどうかで判定）
    if logger.isEnabledFor(DEBUG):
        logger.debug(f"[{category}.{plugin_name}] parse_options結果:")
        logger.debug(f"  入力オプション: {current_options}")
        logger.debug(f"  処理したオプション: {processed}")
        for key, value in processed.items():
            logger.debug(f"    {key}: {value} (型: {type(value).__name__})")
        logger.debug(f"  処理しなかったオプション: {list(unprocessed.keys())}")

    return processed, unprocessed, True


def execute_plugin_chain(
    options_dict: Dict[str, Any],
) -> Tuple[Dict[str, Any], List[Tuple[str, Dict[str, Any], Dict[str, Any]]]]:
    """
    プラグインチェーンを動的に実行

    Args:
        options_dict: 階層的なオプション辞書

    Returns:
        (options_dict, plugin_results) のタプル
        - options_dict: 更新されたオプション辞書（baseオプションから処理済みオプションが削除される）
        - plugin_results: プラグイン実行結果のリスト [(category.name, processed, unprocessed), ...]
    """
    genice3_options = options_dict.get("genice3", {}).copy()
    # debugオプションが指定されていればloggerのlevelをDEBUGにする
    # 指定されていない場合はINFOレベルのまま（既にgenice.pyで設定済み）
    if genice3_options.get("debug"):
        logger.setLevel(DEBUG)
        for handler in logger.handlers:
            handler.setLevel(DEBUG)
    plugin_results: List[Tuple[str, Dict[str, Any], Dict[str, Any]]] = []

    # options_dictをコピーして作業
    updated_options_dict = options_dict.copy()
    updated_options_dict["genice3"] = genice3_options

    for plugin_type in updated_options_dict.keys():
        if plugin_type == "genice3":
            continue
        plugin_config = updated_options_dict[plugin_type]
        if not isinstance(plugin_config, dict):
            continue
        plugin_name = plugin_config.get("name")
        if not plugin_name:
            continue

        # まずgenice3内のオプションを処理し、消費したオプションは消す。
        processed_genice3_options, unprocessed_genice3_options, _ = execute_plugin(
            plugin_type, plugin_name, genice3_options
        )
        plugin_results.append(
            (
                f"{plugin_type}.{plugin_name}",
                processed_genice3_options,
                unprocessed_genice3_options,
            )
        )
        for option in processed_genice3_options:
            genice3_options.pop(option, None)
        logger.debug(
            f"genice3オプションのうち{plugin_type}.{plugin_name}に処理されたオプション: {list(processed_genice3_options.keys())}"
        )
        logger.debug(
            f"genice3オプションのうち{plugin_type}.{plugin_name}に処理されなかったオプション: {list(unprocessed_genice3_options.keys())}"
        )

        # プラグインのオプションを処理（nameを除く）
        plugin_options = {k: v for k, v in plugin_config.items() if k != "name"}
        processed_options, unprocessed_options, _ = execute_plugin(
            plugin_type, plugin_name, plugin_options
        )
        plugin_results.append(
            (f"{plugin_type}.{plugin_name}", processed_options, unprocessed_options)
        )
        logger.debug(
            f"プラグイン{plugin_type}.{plugin_name}のオプション: {plugin_options}"
        )
        logger.debug(
            f"プラグインオプションのうち{plugin_type}.{plugin_name}に処理されたオプション: {list(processed_options.keys())}"
        )
        logger.debug(
            f"プラグインオプションのうち{plugin_type}.{plugin_name}に処理されなかったオプション: {list(unprocessed_options.keys())}"
        )

    updated_options_dict["genice3"] = genice3_options
    return updated_options_dict, plugin_results


# ============================================================================
# メインのパース関数
# ============================================================================


def parse_args(
    args: List[str],
) -> Tuple[Dict[str, Any], List[Tuple[str, Dict[str, Any], Dict[str, Any]]]]:
    """
    コマンドライン引数をパース

    Args:
        args: コマンドライン引数のリスト（sys.argv[1:]相当）

    Returns:
        (options_dict, plugin_results) のタプル
        - options_dict: 階層的なオプション辞書
        - plugin_results: プラグイン実行結果のリスト

    注意: --configオプションが指定されている場合は、最初に処理される必要があります。
    """
    # まず--configオプションを探して、設定ファイルを先に読み込む
    config_dict: Dict[str, Any] = {}
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "--config" or arg == "-C":
            i += 1
            if i < len(args):
                config_path = args[i]
                # 設定ファイルを読み込んで階層的な辞書として取得
                config_dict = load_config_file(config_path)
            i += 1
        else:
            i += 1

    logger.debug(f"config_dict: {config_dict}")

    # コマンドライン引数をフラットな辞書に変換
    cmdline_dict, cmdline_specified_keys = parse_command_line_to_dict(args)

    # デバッグモードで、コマンドライン引数パース直後の状態を表示
    if logger.isEnabledFor(DEBUG):
        logger.debug("=" * 60)
        logger.debug("コマンドライン引数パース直後の状態:")
        logger.debug(f"  config_dict: {config_dict}")
        logger.debug(f"  cmdline_dict: {cmdline_dict}")
        logger.debug(f"  cmdline_specified_keys: {cmdline_specified_keys}")
        logger.debug("=" * 60)

    # オプション辞書を構築（階層的な辞書を作成）
    options_dict = build_options_dict(config_dict, cmdline_dict, cmdline_specified_keys)

    # プラグインチェーンを動的に実行
    options_dict, plugin_results = execute_plugin_chain(options_dict)

    return options_dict, plugin_results


def build_result(
    options_dict: Dict[str, Any],
    plugin_results: List[Tuple[str, Dict[str, Any], Dict[str, Any]]],
) -> Dict[str, Any]:
    """
    パース結果を構築

    Args:
        options_dict: 階層的なオプション辞書
        plugin_results: プラグイン実行結果のリスト

    Returns:
        パース結果の辞書
    """
    # プラグインチェーンの実行結果から、unitcellとexporterの処理結果を取得
    unitcell_processed = {}
    exporter_processed = {}

    # 最後の処理結果を使用（同じプラグインが複数回実行される場合）
    for plugin_name, processed, unprocessed in plugin_results:
        if plugin_name.startswith("unitcell."):
            unitcell_processed.update(processed)
        elif plugin_name.startswith("exporter."):
            exporter_processed.update(processed)

    result = {
        "base_options": options_dict.get(
            "genice3", {}
        ),  # 後方互換性のためbase_optionsという名前で返す
        "unitcell": {
            "name": options_dict.get("unitcell", {}).get("name"),
            "options": {
                k: v for k, v in options_dict.get("unitcell", {}).items() if k != "name"
            },
            "processed": unitcell_processed,
        },
        "exporter": {
            "name": options_dict.get("exporter", {}).get("name"),
            "options": {
                k: v for k, v in options_dict.get("exporter", {}).items() if k != "name"
            },
            "processed": exporter_processed,
        },
        "plugin_chain": [
            {"name": name, "processed": proc, "unprocessed": unproc}
            for name, proc, unproc in plugin_results
        ],
    }
    return result


def validate_result(options_dict: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """
    パース結果を検証

    Args:
        options_dict: 階層的なオプション辞書

    Returns:
        (is_valid, error_messages)
    """
    errors = []

    # unitcell名が必須
    unitcell_name = options_dict.get("unitcell", {}).get("name")
    if not unitcell_name:
        errors.append("unitcell名が指定されていません")

    return len(errors) == 0, errors


# ============================================================================
# 後方互換性のためのPoolBasedParserクラス
# ============================================================================


class PoolBasedParser:
    """
    プールベースのオプションパーサー（後方互換性のためのラッパークラス）

    このクラスは、新しい関数ベースのAPIのラッパーとして機能します。
    新しいコードでは、直接関数（parse_args, build_result, validate_result）を使用することを推奨します。
    """

    def __init__(self):
        self.options_dict: Dict[str, Any] = {}
        self.plugin_results: List[Tuple[str, Dict[str, Any], Dict[str, Any]]] = []

    def parse_args(self, args: List[str]) -> None:
        """
        コマンドライン引数をパース

        Args:
            args: コマンドライン引数のリスト（sys.argv[1:]相当）

        注意: --configオプションが指定されている場合は、最初に処理される必要があります。
        """
        self.options_dict, self.plugin_results = parse_args(args)

    def get_result(self) -> Dict[str, Any]:
        """
        パース結果を取得

        Returns:
            パース結果の辞書
        """
        return build_result(self.options_dict, self.plugin_results)

    def validate(self) -> Tuple[bool, List[str]]:
        """
        パース結果を検証

        Returns:
            (is_valid, error_messages)
        """
        return validate_result(self.options_dict)
