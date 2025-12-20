"""
設定ファイルとコマンドラインの値の上書きをテストする
"""

import sys
import os
import tempfile

# パスを追加
sys.path.insert(0, os.path.dirname(__file__))

from pool_parser import PoolBasedParser
from A15 import parse_options as parse_a15_options
from gromacs import parse_options as parse_gromacs_options


def test_config_overridden_by_cmdline():
    """設定ファイルの値がコマンドラインの値で上書きされることを確認"""
    print("=" * 60)
    print("設定ファイルとコマンドラインの上書きテスト")
    print("=" * 60)

    # 設定ファイルを作成
    yaml_content = """
unitcell:
  name: A15
  shift: [0.2, 0.2, 0.2]
  density: 0.9
  anion:
    "15": Cl
    "16": F
  cation:
    "21": Na

exporter:
  name: gromacs
  guest:
    A12: me
    A13: pr
  water_model: "tip4p"
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        parser = PoolBasedParser()
        # コマンドラインで一部のオプションを指定
        args = [
            "--config",
            config_path,
            "A15",
            "--shift",
            "0.1",
            "0.1",
            "0.1",
            "--anion",
            "17=Br",  # コマンドラインで指定（設定ファイルの値は上書きされる）
            "--density",
            "0.8",
        ]
        parser.parse_args(args)

        result = parser.get_result()
        print("\nunitcellオプション:")
        for key, value in result["unitcell"]["options"].items():
            print(f"  {key}: {value}")

        # A15プラグインで処理
        a15_processed, a15_unprocessed = parse_a15_options(
            result["unitcell"]["options"]
        )

        print("\nA15プラグインが処理したオプション:")
        for key, value in a15_processed.items():
            print(f"  {key}: {value}")

        # 確認: コマンドラインで指定された値が使用されている
        assert a15_processed["shift"] == [0.1, 0.1, 0.1], "shiftはコマンドラインの値"
        assert a15_processed["density"] == 0.8, "densityはコマンドラインの値"
        # anionはコマンドラインで指定された値のみ（設定ファイルの値は上書きされる）
        assert "17" in a15_processed["anion"], "anionにコマンドラインの値が含まれる"
        assert "15" not in a15_processed["anion"], "設定ファイルのanion値は上書きされる"
        assert "16" not in a15_processed["anion"], "設定ファイルのanion値は上書きされる"
        # cationはコマンドラインで指定されていないので、設定ファイルの値が使用される
        assert "21" in a15_processed["cation"], "cationは設定ファイルの値が使用される"

        print("\n✓ テスト成功: コマンドラインの値が設定ファイルの値を上書きしました")

    finally:
        os.unlink(config_path)


def test_multiple_anion_specifications():
    """複数回指定されたanionオプションのテスト"""
    print("\n" + "=" * 60)
    print("複数回指定されたanionオプションのテスト")
    print("=" * 60)

    # 設定ファイルを作成
    yaml_content = """
unitcell:
  name: A15
  anion:
    "15": Cl
    "16": F
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
        f.write(yaml_content)
        config_path = f.name

    try:
        parser = PoolBasedParser()
        # コマンドラインで複数のanionを指定
        args = [
            "--config",
            config_path,
            "A15",
            "--anion",
            "17=Br",
            "--anion",
            "18=I",  # 複数回指定
        ]
        parser.parse_args(args)

        result = parser.get_result()
        print("\nunitcellオプション:")
        print(f"  anion: {result['unitcell']['options']['anion']}")

        # A15プラグインで処理
        a15_processed, a15_unprocessed = parse_a15_options(
            result["unitcell"]["options"]
        )

        print("\nA15プラグインが処理したオプション:")
        print(f"  anion: {a15_processed['anion']}")

        # 確認: コマンドラインで指定された値のみ（設定ファイルの値は上書きされる）
        assert "17" in a15_processed["anion"], "anionにコマンドラインの値が含まれる"
        assert "18" in a15_processed["anion"], "anionにコマンドラインの値が含まれる"
        assert "15" not in a15_processed["anion"], "設定ファイルのanion値は上書きされる"
        assert "16" not in a15_processed["anion"], "設定ファイルのanion値は上書きされる"

        print("\n✓ テスト成功: コマンドラインで複数指定された値が正しく処理されました")

    finally:
        os.unlink(config_path)


if __name__ == "__main__":
    try:
        import yaml
    except ImportError:
        print("YAMLモジュールが利用できません。テストをスキップします。")
        sys.exit(0)

    test_config_overridden_by_cmdline()
    test_multiple_anion_specifications()
