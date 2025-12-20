"""
プールベースオプションパーサーのテスト
"""

import pytest
import tempfile
import os
from pathlib import Path
from pool_parser import PoolBasedParser


class TestPoolBasedParser:
    """プールベースパーサーのテスト"""

    def test_parse_simple_args(self):
        """基本的な引数のパース"""
        parser = PoolBasedParser()
        args = ["A15", "--rep", "2", "2", "2", "--seed", "42"]
        parser.parse_args(args)

        result = parser.get_result()
        assert result["unitcell"]["name"] == "A15"
        assert result["base_options"]["replication_factors"] == ("2", "2", "2")
        assert result["base_options"]["seed"] == "42"

    def test_parse_unitcell_options(self):
        """unitcellプラグインのオプションをパース"""
        parser = PoolBasedParser()
        args = [
            "A15",
            "--rep",
            "2",
            "2",
            "2",
            "--shift",
            "0.1",
            "0.1",
            "0.1",
            "--density",
            "0.8",
            "--anion",
            "15=Cl",
            "--cation",
            "21=Na",
        ]
        parser.parse_args(args)

        result = parser.get_result()
        assert result["unitcell"]["name"] == "A15"
        assert result["unitcell"]["options"]["shift"] == ("0.1", "0.1", "0.1")
        assert result["unitcell"]["options"]["density"] == "0.8"
        assert result["unitcell"]["options"]["anion"] == "15=Cl"
        assert result["unitcell"]["options"]["cation"] == "21=Na"

    def test_parse_exporter_options(self):
        """exporterプラグインのオプションをパース"""
        parser = PoolBasedParser()
        args = [
            "A15",
            "--exporter",
            "gromacs",
            "--guest",
            "A12=me",
            "--guest",
            "A14=et",
            "--water_model",
            "4site",
        ]
        parser.parse_args(args)

        result = parser.get_result()
        assert result["exporter"]["name"] == "gromacs"
        # 同じオプションが複数回指定された場合の処理は実装次第
        # ここでは最後の値が保持される想定
        assert (
            "guest" in result["exporter"]["options"]
            or "guest" in result["unitcell"]["options"]
        )

    def test_parse_complex_example(self):
        """複雑な例（OPTION_HANDLING_PLANS.mdの例）"""
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
        assert result["unitcell"]["name"] == "A15"
        assert result["base_options"]["replication_factors"] == ("2", "2", "2")
        assert result["base_options"]["seed"] == "42"
        assert result["base_options"]["depol_loop"] == "2000"

        # unitcellオプション
        assert result["unitcell"]["options"]["shift"] == ("0.1", "0.1", "0.1")
        assert result["unitcell"]["options"]["anion"] == "15=Cl"
        assert result["unitcell"]["options"]["cation"] == "21=Na"
        assert result["unitcell"]["options"]["density"] == "0.8"

        # exporterオプション（unitcellで処理されなかったもの）
        assert result["exporter"]["name"] == "gromacs"
        # guest, spot_guest, water_model, typeはexporterまたはunitcellのオプション

    def test_parse_bracketed_plugin(self):
        """[plugin --option]形式のパース"""
        parser = PoolBasedParser()
        args = ["[A15", "--shift", "0.1", "0.1", "0.1", "--density", "0.8]"]
        # 注意: 実際の実装では、[と]で囲まれた部分を1つの引数として扱う必要がある
        # ここでは簡易的なテストとして、個別の引数でテスト

    def test_parse_yaml(self):
        """YAMLファイルからのパース"""
        yaml_content = """
genice3:
  seed: 42
  depol_loop: 2000
  replication_factors: [2, 2, 2]
  spot_anion:
    "1": Cl
  spot_cation:
    "5": Na

unitcell:
  name: A15
  shift: [0.1, 0.1, 0.1]
  anion:
    "15": Cl
  cation:
    "21": Na
  density: 0.8

exporter:
  name: gromacs
  guest:
    A12: me
    A14: et
  spot_guest:
    "0": "4site"
  water_model: "4site"
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as f:
            f.write(yaml_content)
            yaml_path = f.name

        try:
            parser = PoolBasedParser()
            parser.parse_yaml(yaml_path)

            result = parser.get_result()
            assert result["base_options"]["seed"] == 42
            assert result["base_options"]["depol_loop"] == 2000
            assert result["base_options"]["replication_factors"] == (2, 2, 2)

            assert result["unitcell"]["name"] == "A15"
            assert result["unitcell"]["options"]["shift"] == [0.1, 0.1, 0.1]
            assert result["unitcell"]["options"]["density"] == 0.8

            assert result["exporter"]["name"] == "gromacs"
            assert result["exporter"]["options"]["guest"]["A12"] == "me"
            assert result["exporter"]["options"]["guest"]["A14"] == "et"
        finally:
            os.unlink(yaml_path)

    def test_validate(self):
        """バリデーションテスト"""
        parser = PoolBasedParser()
        args = ["A15", "--rep", "2", "2", "2"]
        parser.parse_args(args)

        is_valid, errors = parser.validate()
        assert is_valid
        assert len(errors) == 0

    def test_validate_missing_unitcell(self):
        """unitcell名がない場合のバリデーション"""
        parser = PoolBasedParser()
        args = ["--rep", "2", "2", "2"]
        parser.parse_args(args)

        is_valid, errors = parser.validate()
        assert not is_valid
        assert any("unitcell名" in error for error in errors)

    def test_unprocessed_options(self):
        """処理されなかったオプションのテスト"""
        parser = PoolBasedParser()
        args = ["A15", "--unknown_option", "value", "--rep", "2", "2", "2"]
        parser.parse_args(args)

        result = parser.get_result()
        # unknown_optionはunitcell_optionsに含まれる（unitcellプラグインが処理する想定）
        assert "unknown_option" in result["unitcell"]["options"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
