#!/usr/bin/env python3
"""
コマンドライン引数をパースする簡易スクリプト

使用方法:
    python parse_args.py A15 --exporter gromacs --water_model foursite --type ice
    python parse_args.py A15 --rep 2 2 2 --shift 0.1 0.1 0.1 --density 0.8
    python parse_args.py --config test_config.yaml A15 --exporter gromacs
"""

import sys
from pool_parser import PoolBasedParser
import json


def format_value(value):
    """値を読みやすい形式でフォーマット"""
    if isinstance(value, dict):
        return json.dumps(value, indent=2, ensure_ascii=False)
    elif isinstance(value, (list, tuple)):
        return str(list(value))
    else:
        return str(value)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    # コマンドライン引数を取得（スクリプト名を除く）
    args = sys.argv[1:]

    # パーサーを作成して実行
    parser = PoolBasedParser()
    try:
        parser.parse_args(args)
        result = parser.get_result()

        # 結果を表示
        print("=" * 60)
        print("パース結果")
        print("=" * 60)

        print("\n基底レベルのオプション:")
        if result["base_options"]:
            for key, value in result["base_options"].items():
                print(f"  {key}: {format_value(value)}")
        else:
            print("  （なし）")

        print(f"\nunitcell: {result['unitcell']['name'] or '(指定なし)'}")
        if result["unitcell"]["options"]:
            print("unitcellオプション（入力）:")
            for key, value in result["unitcell"]["options"].items():
                print(f"  {key}: {format_value(value)}")
        else:
            print("  （なし）")

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
                        print(f"    {key}: {format_value(value)}")
                else:
                    print("  処理したオプション: （なし）")
                if unprocessed:
                    print("  処理しなかったオプション:")
                    for key, value in unprocessed.items():
                        print(f"    {key}: {format_value(value)}")
                else:
                    print("  処理しなかったオプション: （なし）")

        print(f"\nexporter: {result['exporter']['name'] or '(指定なし)'}")
        if result["exporter"]["options"]:
            print("exporterオプション（設定ファイルからの値）:")
            for key, value in result["exporter"]["options"].items():
                print(f"  {key}: {format_value(value)}")
        else:
            print("  （なし）")

        if result["unprocessed_options"]:
            print("\n処理されなかったオプション:")
            for key, value in result["unprocessed_options"].items():
                print(f"  {key}: {format_value(value)}")
        else:
            print("\n処理されなかったオプション: （なし）")

        # バリデーション
        is_valid, errors = parser.validate()
        if is_valid:
            print("\n✓ バリデーション成功")
        else:
            print("\n✗ バリデーションエラー:")
            for error in errors:
                print(f"  - {error}")
            sys.exit(1)

    except Exception as e:
        print(f"エラー: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
