"""
プールベースオプションパーサーの使用例
"""

from pool_parser import PoolBasedParser


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
        "0=4site",
        "--water_model",
        "4site",
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
