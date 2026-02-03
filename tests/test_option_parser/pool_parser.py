"""
プールベースのオプションパーサー

基底レベルのオプションを処理し、処理されなかったオプションを
プラグインに順次渡していく方式のパーサー。
"""

from typing import Dict, Any, List, Optional, Set, Tuple
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
        self.exporter_name: Optional[str] = None
        self.exporter_options: Dict[str, Any] = {}
        self.unprocessed_options: Dict[str, Any] = {}

    def parse_args(self, args: List[str]) -> None:
        """
        コマンドライン引数をパース

        Args:
            args: コマンドライン引数のリスト（sys.argv[1:]相当）
        """
        pool = OptionPool({})
        i = 0
        unitcell_arg = None

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

            # --で始まるオプション
            if arg.startswith("--"):
                key = arg[2:]  # --を削除
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
            self.unitcell_options = unprocessed.copy()
            # unitcellプラグインが処理しないオプションはexporterに渡す想定
            # ここでは、unitcellオプションとして全て保持しておく
            # （実際のプラグイン処理で、処理されなかったものがexporterに渡される）
        else:
            # unitcell名が指定されていない場合は、すべての未処理オプションを保持
            self.unprocessed_options = unprocessed.copy()

        # 注意: 実際の実装では、unitcellプラグインがオプションを処理した後、
        # 処理されなかったオプションをexporterプラグインに渡す必要がある

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
        unprocessed = pool.get_unprocessed()
        if self.unitcell_name:
            self.unitcell_options = unprocessed.copy()

        # exporter名を設定
        if self.exporter_name is None and "exporter" in pool:
            self.exporter_name = pool.get("exporter")

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
        result = {
            "base_options": self.base_options,
            "unitcell": {
                "name": self.unitcell_name,
                "options": self.unitcell_options,
            },
            "exporter": {
                "name": self.exporter_name,
                "options": self.exporter_options,
            },
            "unprocessed_options": self.unprocessed_options,
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
