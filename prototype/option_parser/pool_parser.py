"""
プールベースのオプションパーサー

基底レベルのオプションを処理し、処理されなかったオプションを
プラグインに順次渡していく方式のパーサー。
"""

from typing import Dict, Any, List, Optional, Set, Tuple, Callable
from pathlib import Path

# YAMLライブラリのインポート
try:
    import yaml

    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False


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


class PoolBasedParser:
    """
    プールベースのオプションパーサー

    基底レベルのオプションを処理し、処理されなかったオプションを
    プラグインに順次渡していく。
    """

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

    def __init__(self):
        self.base_options: Dict[str, Any] = {}
        self.unitcell_name: Optional[str] = None
        self.unitcell_options: Dict[str, Any] = {}
        self.unitcell_config_options: Dict[str, Any] = (
            {}
        )  # 設定ファイルからのunitcellオプション
        self.exporter_name: Optional[str] = None
        self.exporter_options: Dict[str, Any] = {}
        self.exporter_config_options: Dict[str, Any] = (
            {}
        )  # 設定ファイルからのexporterオプション
        self.unprocessed_options: Dict[str, Any] = {}
        # 動的に実行されたプラグインの処理結果
        self.plugin_results: List[Tuple[str, Dict[str, Any], Dict[str, Any]]] = []
        # プラグイン名 -> (category, parse_function) のマッピング
        self._plugin_cache: Dict[Tuple[str, str], Callable] = {}

    def parse_args(self, args: List[str]) -> None:
        """
        コマンドライン引数をパース

        Args:
            args: コマンドライン引数のリスト（sys.argv[1:]相当）

        注意: --configオプションが指定されている場合は、最初に処理される必要があります。
        """
        pool = OptionPool({})
        i = 0
        unitcell_arg = None
        config_path = None

        # まず--configオプションを探して、設定ファイルを先に読み込む
        while i < len(args):
            arg = args[i]
            if arg == "--config" or arg == "-C":
                i += 1
                if i < len(args):
                    config_path = args[i]
                    # 設定ファイルを読み込む
                    self._load_config_file(config_path)
                i += 1
            else:
                i += 1

        # 再度最初からループして、通常のオプションを処理
        # コマンドラインで指定されたキーを追跡（設定ファイルの値を上書きするため）
        cmdline_specified_keys: Set[str] = set()
        i = 0
        while i < len(args):
            arg = args[i]

            # unitcell名の位置引数を処理
            if not arg.startswith("--") and unitcell_arg is None:
                unitcell_arg = arg
                # [plugin --option]形式をチェック
                if arg.startswith("[") and arg.endswith("]"):
                    unitcell_name, options = self._parse_bracketed_plugin(arg)
                    self.unitcell_name = unitcell_name
                    self._add_to_pool(pool, options)
                elif "[" in arg and "]" in arg:
                    # plugin[--option value]形式
                    left = arg.find("[")
                    right = arg.rfind("]")
                    unitcell_name = arg[:left]
                    option_str = arg[left + 1 : right]
                    self.unitcell_name = unitcell_name
                    self._add_to_pool(pool, self._parse_option_string(option_str))
                else:
                    self.unitcell_name = arg
                i += 1
                continue

            # --configオプションは既に処理済みなのでスキップ
            if arg == "--config" or arg == "-C":
                i += 2  # --config とその値をスキップ
                continue

            # --で始まるオプション
            if arg.startswith("--"):
                key = arg[2:]  # --を削除
                # コマンドラインで指定されたキーとして記録
                cmdline_specified_keys.add(key)
                i += 1

                # 値を取得（次の引数が--で始まっていなければ値として扱う）
                values: List[Any] = []
                while i < len(args) and not args[i].startswith("--"):
                    values.append(args[i])
                    i += 1

                # 値が1つの場合はそのまま、複数の場合はタプル
                value = values[0] if len(values) == 1 else tuple(values)

                # 同じオプションが複数回指定された場合はリストに変換
                if key in pool.options:
                    existing = pool.options[key]
                    if not isinstance(existing, list):
                        pool.options[key] = [existing]
                    pool.options[key].append(value)
                else:
                    pool.options[key] = value
            else:
                i += 1

        # 基底レベルのオプションを処理
        self._parse_base_level(pool)

        # unitcell名がまだ設定されていない場合、poolから取得
        if self.unitcell_name is None and "unitcell" in pool:
            self.unitcell_name = pool.get("unitcell")

        # exporter名を設定
        if "exporter" in pool:
            self.exporter_name = pool.get("exporter")

        # 残りのオプションをunitcellプラグインに渡す
        unprocessed = pool.get_unprocessed()
        if self.unitcell_name:
            # 設定ファイルのオプションを初期値として使用
            self.unitcell_options = self.unitcell_config_options.copy()
            # コマンドラインオプションで上書き
            # コマンドラインで指定されたキーは設定ファイルの値を完全に上書き
            for key, value in unprocessed.items():
                if key in cmdline_specified_keys:
                    # コマンドラインで指定された場合は、設定ファイルの値を完全に上書き
                    self.unitcell_options[key] = value
                # コマンドラインで指定されていない場合は、設定ファイルの値を使用（既に設定済み）

        # プラグインチェーンを動的に実行
        self._execute_plugin_chain()

    def _load_plugin(self, category: str, name: str) -> Optional[Callable]:
        """
        プラグインを動的に読み込む

        Args:
            category: プラグインのカテゴリ（"unitcell", "exporter", "molecule"）
            name: プラグイン名

        Returns:
            parse_options関数、またはNone（読み込めない場合）
        """
        import importlib
        import os

        cache_key = (category, name)
        if cache_key in self._plugin_cache:
            return self._plugin_cache[cache_key]

        try:
            # prototype/option_parserディレクトリから直接インポート
            # 現在のディレクトリを取得
            current_dir = os.path.dirname(os.path.abspath(__file__))
            # モジュール名は単純にname（ファイル名から拡張子を除いたもの）
            module = importlib.import_module(name)
            if hasattr(module, "parse_options"):
                parse_func = module.parse_options
                self._plugin_cache[cache_key] = parse_func
                return parse_func
        except (ImportError, AttributeError, ModuleNotFoundError):
            pass

        return None

    def _execute_plugin_chain(self) -> None:
        """
        プラグインチェーンを動的に実行
        unitcell → exporter → molecule (water_modelに基づく) の順で実行
        """
        current_options = self.unitcell_options.copy()
        self.plugin_results = []

        # 1. unitcellプラグインを実行
        if self.unitcell_name:
            unitcell_parse = self._load_plugin("unitcell", self.unitcell_name)
            if unitcell_parse:
                processed, unprocessed = unitcell_parse(current_options)
                self.plugin_results.append(
                    (f"unitcell.{self.unitcell_name}", processed, unprocessed)
                )
                # exporterオプションと統合
                current_options = {**self.exporter_config_options, **unprocessed}

        # 2. exporterプラグインを実行
        if self.exporter_name:
            exporter_parse = self._load_plugin("exporter", self.exporter_name)
            if exporter_parse:
                processed, unprocessed = exporter_parse(current_options)
                self.plugin_results.append(
                    (f"exporter.{self.exporter_name}", processed, unprocessed)
                )

                # water_modelが処理された場合、その値に基づいてmoleculeプラグインをチェーンに追加
                if "water_model" in processed:
                    water_model = processed["water_model"]
                    # water_modelの値を取得（文字列または辞書から）
                    water_model_name = (
                        water_model.get("name")
                        if isinstance(water_model, dict)
                        else water_model
                    )
                    if isinstance(water_model_name, str):
                        # moleculeプラグインを動的にチェーンに追加
                        molecule_parse = self._load_plugin("molecule", water_model_name)
                        if molecule_parse:
                            # 処理済みのオプションも含めてmoleculeに渡す
                            # （water_modelの値を確認するために必要）
                            molecule_options = {**processed, **unprocessed}
                            molecule_processed, molecule_unprocessed = molecule_parse(
                                molecule_options
                            )
                            self.plugin_results.append(
                                (
                                    f"molecule.{water_model_name}",
                                    molecule_processed,
                                    molecule_unprocessed,
                                )
                            )
                            # moleculeで処理したオプションをexporterの処理結果に統合
                            # 最後に追加した結果を更新
                            last_name, last_processed, last_unprocessed = (
                                self.plugin_results[-2]
                            )  # exporterの結果
                            if molecule_processed:
                                if isinstance(last_processed["water_model"], str):
                                    last_processed["water_model"] = {
                                        "name": last_processed["water_model"],
                                        "options": molecule_processed,
                                    }
                                elif isinstance(last_processed["water_model"], dict):
                                    if "options" not in last_processed["water_model"]:
                                        last_processed["water_model"]["options"] = {}
                                    last_processed["water_model"]["options"].update(
                                        molecule_processed
                                    )
                            self.plugin_results[-2] = (
                                last_name,
                                last_processed,
                                last_unprocessed,
                            )
                            # moleculeで処理したオプションはunprocessedから除外
                            # （moleculeで処理されたキーをunprocessedから削除）
                            current_options = {
                                k: v
                                for k, v in unprocessed.items()
                                if k not in molecule_processed
                            }
                        else:
                            current_options = unprocessed
                    else:
                        current_options = unprocessed
                else:
                    current_options = unprocessed
            else:
                current_options = unprocessed

        # 最終的に処理されなかったオプション（どのプラグインでも処理されなかったもの）
        self.unprocessed_options = current_options

    def _load_config_file(self, yaml_path: str) -> None:
        """
        設定ファイルを読み込む（内部メソッド）

        Args:
            yaml_path: YAMLファイルのパス
        """
        if not YAML_AVAILABLE:
            return  # YAMLが利用できない場合は何もしない

        config_file = Path(yaml_path)
        if not config_file.exists():
            return  # ファイルが存在しない場合は何もしない

        with open(config_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}

        # unitcellセクション
        if "unitcell" in config:
            unitcell_config = config["unitcell"]
            if isinstance(unitcell_config, dict):
                if "name" in unitcell_config:
                    # コマンドラインで指定されていない場合のみ設定ファイルの名前を使用
                    if self.unitcell_name is None:
                        self.unitcell_name = unitcell_config["name"]
                # unitcellのオプションを保存（nameを除く）
                self.unitcell_config_options = {
                    k: v for k, v in unitcell_config.items() if k != "name"
                }

        # exporterセクション
        if "exporter" in config:
            exporter_config = config["exporter"]
            if isinstance(exporter_config, dict):
                if "name" in exporter_config:
                    # コマンドラインで指定されていない場合のみ設定ファイルの名前を使用
                    if self.exporter_name is None:
                        self.exporter_name = exporter_config["name"]
                # exporterのオプションを保存（nameを除く）
                self.exporter_config_options = {
                    k: v for k, v in exporter_config.items() if k != "name"
                }

        # genice3セクション（基底レベルのオプション）も処理
        if "genice3" in config:
            genice3_config = config["genice3"]
            # 基底レベルのオプションは、後でコマンドライン引数で上書きされる

    def parse_yaml(self, yaml_path: str) -> None:
        """
        YAML設定ファイルをパース

        Args:
            yaml_path: YAMLファイルのパス
        """
        if not YAML_AVAILABLE:
            raise ValueError(
                "YAMLライブラリが利用できません。`pip install pyyaml`でインストールしてください。"
            )

        config_file = Path(yaml_path)
        if not config_file.exists():
            raise FileNotFoundError(f"設定ファイルが見つかりません: {yaml_path}")

        with open(config_file, "r", encoding="utf-8") as f:
            config = yaml.safe_load(f) or {}

        pool = OptionPool({})

        # genice3セクション（基底レベルのオプション）
        if "genice3" in config:
            genice3_config = config["genice3"]
            for key, value in genice3_config.items():
                pool.options[key] = value

        # unitcellセクション
        if "unitcell" in config:
            unitcell_config = config["unitcell"]
            if isinstance(unitcell_config, dict):
                if "name" in unitcell_config:
                    self.unitcell_name = unitcell_config["name"]
                # unitcellのオプションをプールに追加
                for key, value in unitcell_config.items():
                    if key != "name":
                        pool.options[key] = value

        # exporterセクション
        if "exporter" in config:
            exporter_config = config["exporter"]
            if isinstance(exporter_config, dict):
                if "name" in exporter_config:
                    self.exporter_name = exporter_config["name"]
                    pool.options["exporter"] = exporter_config["name"]
                # exporterのオプションをプールに追加
                for key, value in exporter_config.items():
                    if key != "name":
                        pool.options[key] = value

        # 基底レベルのオプションを処理
        self._parse_base_level(pool)

        # 残りのオプションをunitcellプラグインに渡す
        # 設定ファイルのオプションを初期値として、コマンドラインオプションで上書き
        unprocessed = pool.get_unprocessed()
        if self.unitcell_name:
            # 設定ファイルのオプションを初期値として使用
            self.unitcell_options = self.unitcell_config_options.copy()
            # コマンドラインオプションで上書き
            self.unitcell_options.update(unprocessed)

        # exporter名を設定
        if self.exporter_name is None and "exporter" in pool:
            self.exporter_name = pool.get("exporter")

        # exporterオプションも設定ファイルの値を初期値として使用
        # （ただし、parse_yamlではコマンドライン引数がないので、設定ファイルの値のみ）
        if self.exporter_name:
            self.exporter_options = self.exporter_config_options.copy()

    def _parse_base_level(self, pool: OptionPool) -> None:
        """基底レベルのオプションを処理"""
        for key in self.BASE_LEVEL_OPTIONS:
            if key in pool:
                value = pool.get(key)
                if key == "rep":
                    # repはreplication_factorsのエイリアス
                    self.base_options["replication_factors"] = (
                        value if isinstance(value, tuple) else tuple(value)
                    )
                elif key == "spot_anion" or key == "spot_cation":
                    # key=value形式の文字列を辞書に変換
                    if isinstance(value, list):
                        # 複数回指定された場合
                        result_dict = {}
                        for item in value:
                            if isinstance(item, str) and "=" in item:
                                k, v = item.split("=", 1)
                                result_dict[k] = v
                            else:
                                result_dict[str(len(result_dict))] = item
                        self.base_options[key] = result_dict
                    elif isinstance(value, str):
                        if "=" in value:
                            k, v = value.split("=", 1)
                            self.base_options.setdefault(key, {})[k] = v
                        else:
                            self.base_options[key] = value
                    elif isinstance(value, dict):
                        self.base_options[key] = value
                    else:
                        self.base_options[key] = value
                else:
                    self.base_options[key] = value

    def _parse_bracketed_plugin(self, arg: str) -> Tuple[str, Dict[str, Any]]:
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
            options = self._parse_option_string(options_str)
        else:
            options = {}

        return plugin_name, options

    def _parse_option_string(self, option_str: str) -> Dict[str, Any]:
        """
        --option value --option2 value2形式の文字列をパース

        Args:
            option_str: オプション文字列

        Returns:
            オプションの辞書
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
                value = values[0] if len(values) == 1 else tuple(values)
                options[key] = value
            else:
                i += 1
        return options

    def _add_to_pool(self, pool: OptionPool, options: Dict[str, Any]) -> None:
        """オプションをプールに追加"""
        for key, value in options.items():
            pool.options[key] = value

    def get_result(self) -> Dict[str, Any]:
        """
        パース結果を取得

        Returns:
            パース結果の辞書
        """
        # プラグインチェーンの実行結果から、unitcellとexporterの処理結果を取得
        unitcell_processed = {}
        exporter_processed = {}

        for plugin_name, processed, unprocessed in self.plugin_results:
            if plugin_name.startswith("unitcell."):
                unitcell_processed = processed
            elif plugin_name.startswith("exporter."):
                exporter_processed = processed

        result = {
            "base_options": self.base_options,
            "unitcell": {
                "name": self.unitcell_name,
                "options": self.unitcell_options,  # 入力オプション
                "processed": unitcell_processed,  # 処理結果
            },
            "exporter": {
                "name": self.exporter_name,
                "options": self.exporter_config_options.copy(),  # 設定ファイルからのオプション
                "processed": exporter_processed,  # 処理結果
            },
            "unprocessed_options": self.unprocessed_options,
            "plugin_chain": [
                {"name": name, "processed": proc, "unprocessed": unproc}
                for name, proc, unproc in self.plugin_results
            ],
        }
        return result

    def validate(self) -> Tuple[bool, List[str]]:
        """
        パース結果を検証

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # unitcell名が必須
        if not self.unitcell_name:
            errors.append("unitcell名が指定されていません")

        # 処理されていないオプションがある場合は警告（エラーにはしない）
        if self.unprocessed_options:
            errors.append(
                f"処理されなかったオプションがあります: {list(self.unprocessed_options.keys())}"
            )

        return len(errors) == 0, errors
