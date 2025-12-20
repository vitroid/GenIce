"""
プールベースオプションパーサーの使用例
"""

from pool_parser import PoolBasedParser

# プラグインをインポート（利用可能な場合）
try:
    from A15 import parse_options as parse_a15_options
    from gromacs import parse_options as parse_gromacs_options
    from foursite import parse_options as parse_foursite_options

    PLUGINS_AVAILABLE = True
    FOURSITE_AVAILABLE = True
except ImportError:
    PLUGINS_AVAILABLE = False
    FOURSITE_AVAILABLE = False


def example_command_line():
    """コマンドライン引数からのパース例"""
    print("=" * 60)
    print("コマンドライン引数からのパース例")
    print("=" * 60)

    parser = PoolBasedParser()
    args = [
        "A15",
        "--rep",
        "2",
        "2",
        "2",
        "--exporter",
        "gromacs",
        "--seed",
        "42",
        "--depol_loop",
        "2000",
        "--spot_anion",
        "1=Cl",
        "--spot_cation",
        "5=Na",
        "--shift",
        "0.1",
        "0.1",
        "0.1",
        "--anion",
        "15=Cl",
        "--cation",
        "21=Na",
        "--density",
        "0.8",
        "--guest",
        "A12=me",
        "--guest",
        "A14=et",
        "--spot_guest",
        "0=foursite",
        "--water_model",
        "foursite",
        "--type",
        "ice",
    ]
    parser.parse_args(args)

    result = parser.get_result()

    print("\n基底レベルのオプション:")
    for key, value in result["base_options"].items():
        print(f"  {key}: {value}")

    print(f"\nunitcell: {result['unitcell']['name']}")
    print("unitcellオプション:")
    for key, value in result["unitcell"]["options"].items():
        print(f"  {key}: {value}")

    print(f"\nexporter: {result['exporter']['name']}")
    print("exporterオプション（設定ファイルからの値）:")
    for key, value in result["exporter"]["options"].items():
        print(f"  {key}: {value}")

    # 動的に実行されたプラグインチェーンの結果を表示
    if result.get("plugin_chain"):
        print("\n--- 動的に実行されたプラグインチェーン ---")
        for item in result["plugin_chain"]:
            plugin_name = item["name"]
            processed = item["processed"]
            unprocessed = item["unprocessed"]
            print(f"\n{plugin_name}:")
            if processed:
                print("  処理したオプション:")
                for key, value in processed.items():
                    print(f"    {key}: {value}")
            if unprocessed:
                print("  処理しなかったオプション:")
                for key, value in unprocessed.items():
                    print(f"    {key}: {value}")

    # バリデーション
    is_valid, errors = parser.validate()
    if is_valid:
        print("\n✓ バリデーション成功")
    else:
        print("\n✗ バリデーションエラー:")
        for error in errors:
            print(f"  - {error}")


def example_yaml():
    """YAMLファイルからのパース例"""
    try:
        import yaml
    except ImportError:
        print("YAMLモジュールが利用できません。スキップします。")
        return

    print("\n" + "=" * 60)
    print("YAMLファイルからのパース例")
    print("=" * 60)

    parser = PoolBasedParser()
    try:
        parser.parse_yaml("test_config.yaml")

        result = parser.get_result()

        print("\n基底レベルのオプション:")
        for key, value in result["base_options"].items():
            print(f"  {key}: {value}")

        print(f"\nunitcell: {result['unitcell']['name']}")
        print("unitcellオプション:")
        for key, value in result["unitcell"]["options"].items():
            print(f"  {key}: {value}")

        print(f"\nexporter: {result['exporter']['name']}")
        print("exporterオプション:")
        for key, value in result["exporter"]["options"].items():
            print(f"  {key}: {value}")

        # バリデーション
        is_valid, errors = parser.validate()
        if is_valid:
            print("\n✓ バリデーション成功")
        else:
            print("\n✗ バリデーションエラー:")
            for error in errors:
                print(f"  - {error}")
    except FileNotFoundError:
        print("test_config.yamlが見つかりません。スキップします。")


if __name__ == "__main__":
    example_command_line()
    example_yaml()
