#!/usr/bin/env python3
# coding: utf-8
"""
genice2.latticesモジュールをgenice3.unitcellに変換するスクリプト

Usage:
    python convert_to_genice3.py <lattice_file.py> [output_dir]
    python convert_to_genice3.py --all [output_dir]
"""

import ast
import sys
import os
import re
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import importlib.util


def parse_pairs_string(pairs_str: str) -> List[Tuple[int, int]]:
    """pairs文字列をパースしてタプルのリストに変換"""
    pairs = []
    for line in pairs_str.split("\n"):
        cols = line.split()
        if len(cols) == 2:
            try:
                pairs.append((int(cols[0]), int(cols[1])))
            except ValueError:
                pass
    return pairs


def extract_string_value(node: ast.AST) -> Optional[str]:
    """ASTノードから文字列値を抽出"""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    elif hasattr(ast, "Str") and isinstance(node, ast.Str):  # Python < 3.8
        return node.s
    return None


def extract_numeric_value(node: ast.AST) -> Optional[float]:
    """ASTノードから数値を抽出"""
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return float(node.value)
    elif hasattr(ast, "Num") and isinstance(node, ast.Num):  # Python < 3.8
        return float(node.n)
    elif isinstance(node, ast.BinOp):
        # 簡単な計算式（例: 7.78 / 10.0）を評価
        try:
            if hasattr(ast, "unparse"):
                return eval(ast.unparse(node))
            else:
                # Python < 3.9 fallback
                code = compile(ast.Expression(node), "<string>", "eval")
                return eval(code)
        except:
            pass
    return None


def analyze_lattice_file(filepath: Path) -> Dict[str, Any]:
    """latticeファイルを解析して情報を抽出"""
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    try:
        tree = ast.parse(content)
    except SyntaxError as e:
        print(f"Error parsing {filepath}: {e}", file=sys.stderr)
        return {}

    info = {
        "filepath": filepath,
        "module_name": filepath.stem,
        "header": [],
        "imports": [],
        "docstring": None,
        "desc": None,
        "has_cif": False,
        "pairs": None,
        "pairs_str": None,
        "waters": None,
        "waters_str": None,
        "waters_array": None,
        "cell": None,
        "cell_str": None,
        "coord": "relative",
        "bondlen": 3,
        "density": None,
        "fixed": None,
        "fixed_str": None,
        "cages": None,
        "original_content": content,  # 旧コード全体を保存
    }

    # ヘッダー部分（クラス定義より前）を抽出
    lines = content.split("\n")
    in_class = False
    header_end = 0

    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            # クラス定義の行番号を取得
            header_end = node.lineno - 1
            break

    info["header"] = lines[:header_end]

    # docstringを抽出（モジュールレベルのdocstring）
    # ASTから抽出する方法と、正規表現で抽出する方法の両方を試す
    for node in ast.walk(tree):
        if isinstance(node, ast.Module):
            if node.body and isinstance(node.body[0], ast.Expr):
                docstring_value = extract_string_value(node.body[0].value)
                if docstring_value:
                    info["docstring"] = f'"""{docstring_value}"""'
            break

    # 正規表現でも試す（より確実）
    if not info.get("docstring"):
        docstring_match = re.search(r'^"""(.*?)"""', content, re.DOTALL | re.MULTILINE)
        if docstring_match:
            info["docstring"] = f'"""{docstring_match.group(1)}"""'

    # クラス内の__init__メソッドを探す
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                    # __init__メソッド内の属性割り当てを解析
                    for stmt in item.body:
                        if isinstance(stmt, ast.Assign):
                            for target in stmt.targets:
                                if isinstance(target, ast.Attribute):
                                    attr_name = target.attr

                                    # self.pairs
                                    if attr_name == "pairs":
                                        if isinstance(
                                            stmt.value,
                                            (
                                                ast.Constant,
                                                getattr(ast, "Str", type(None)),
                                            ),
                                        ):
                                            info["pairs_str"] = extract_string_value(
                                                stmt.value
                                            )
                                        elif isinstance(stmt.value, ast.Name):
                                            # self.pairs = self.fixed のような場合
                                            info["pairs"] = "fixed"

                                    # self.waters
                                    elif attr_name == "waters":
                                        # 文字列の場合
                                        waters_str = extract_string_value(stmt.value)
                                        if waters_str:
                                            info["waters_str"] = waters_str
                                        # np.array([...])の場合
                                        elif isinstance(stmt.value, ast.Call):
                                            # np.array()呼び出しを検出
                                            if (
                                                isinstance(
                                                    stmt.value.func, ast.Attribute
                                                )
                                                and stmt.value.func.attr == "array"
                                            ) or (
                                                isinstance(stmt.value.func, ast.Name)
                                                and stmt.value.func.id == "array"
                                            ):
                                                # 引数がリストの場合
                                                if stmt.value.args and isinstance(
                                                    stmt.value.args[0], ast.List
                                                ):
                                                    # リストを文字列形式に変換
                                                    if hasattr(ast, "unparse"):
                                                        # ASTを文字列に変換してから評価
                                                        list_str = ast.unparse(
                                                            stmt.value.args[0]
                                                        )
                                                        info["waters_array"] = list_str
                                                    else:
                                                        # Python < 3.9: 手動で抽出
                                                        info["waters_array"] = (
                                                            stmt.value.args[0]
                                                        )
                                                else:
                                                    # その他の場合は検出できず
                                                    pass

                                    # self.cell
                                    elif attr_name == "cell":
                                        if isinstance(stmt.value, ast.Call):
                                            # cellvectors()呼び出し
                                            info["cell"] = (
                                                ast.unparse(stmt.value)
                                                if hasattr(ast, "unparse")
                                                else None
                                            )
                                        else:
                                            info["cell_str"] = extract_string_value(
                                                stmt.value
                                            )

                                    # self.coord
                                    elif attr_name == "coord":
                                        val = extract_string_value(stmt.value)
                                        if val:
                                            info["coord"] = val

                                    # self.bondlen
                                    elif attr_name == "bondlen":
                                        val = extract_numeric_value(stmt.value)
                                        if val is not None:
                                            info["bondlen"] = val

                                    # self.density
                                    elif attr_name == "density":
                                        val = extract_numeric_value(stmt.value)
                                        if val is not None:
                                            info["density"] = val

                                    # self.fixed
                                    elif attr_name == "fixed":
                                        if isinstance(stmt.value, ast.Name):
                                            # self.fixed = self.pairs のような場合
                                            info["fixed"] = "pairs"
                                        else:
                                            info["fixed_str"] = extract_string_value(
                                                stmt.value
                                            )

                                    # self.cages
                                    elif attr_name == "cages":
                                        info["cages"] = extract_string_value(stmt.value)

    # CIF使用の検出
    for line in content.split("\n"):
        if "CIF.waters_and_pairs" in line or "from genice2 import CIF" in line:
            info["has_cif"] = True
            break

    # importsを抽出
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            if hasattr(ast, "unparse"):
                info["imports"].append(ast.unparse(node))
            else:
                # Python < 3.9 fallback
                if isinstance(node, ast.Import):
                    names = ", ".join([alias.name for alias in node.names])
                    info["imports"].append(f"import {names}")
                elif isinstance(node, ast.ImportFrom):
                    names = (
                        ", ".join([alias.name for alias in node.names])
                        if node.names
                        else "*"
                    )
                    info["imports"].append(f"from {node.module} import {names}")

    # desc辞書を抽出（文字列として保持）
    # より正確なマッチング: ネストされた辞書にも対応
    desc_match = re.search(
        r"desc\s*=\s*(\{[^}]*\{[^}]*\}[^}]*\}|\{[^}]+\})", content, re.DOTALL
    )
    if desc_match:
        info["desc"] = desc_match.group(1).strip()
    else:
        # シンプルなパターンも試す
        desc_match = re.search(r"desc\s*=\s*\{[^}]+\}", content, re.DOTALL)
        if desc_match:
            info["desc"] = desc_match.group().replace("desc =", "").strip()

    return info


def generate_unitcell_code(
    info: Dict[str, Any], generate_skeleton: bool = True
) -> Optional[str]:
    """解析結果からgenice3.unitcellコードを生成

    自動変換できない場合でも、generate_skeletonがTrueの場合はスケルトンコードを生成
    """
    lines = []

    # 変換可能かどうかを判定
    can_convert = True
    reason = []

    if info.get("has_cif"):
        can_convert = False
        reason.append("CIFパターンを使用しているため")

    if (
        not info.get("waters_str")
        and not info.get("waters_array")
        and not info.get("has_cif")
    ):
        can_convert = False
        reason.append("watersが定義されていないため")

    if not info.get("cell") and not info.get("cell_str"):
        can_convert = False
        reason.append("cellが定義されていないため")

    # 変換できない場合でもスケルトンを生成
    if not can_convert and not generate_skeleton:
        return None

    # 変換できない場合のスケルトンコード生成
    if not can_convert:
        # 旧コードをコメントとして含める
        if info.get("original_content"):
            lines.append(
                "# ============================================================================"
            )
            lines.append("# Original code from genice2.lattices (for reference)")
            lines.append(
                "# ============================================================================"
            )
            original_lines = info["original_content"].split("\n")
            for line in original_lines:
                lines.append(f"# {line}")
            lines.append("")
            lines.append(
                "# ============================================================================"
            )
            lines.append("# End of original code")
            lines.append(
                "# ============================================================================"
            )
            lines.append("")

        # docstringをそのまま移植
        if info.get("docstring"):
            lines.append(info["docstring"])
            lines.append("")

        # descをそのまま移植
        if info.get("desc"):
            lines.append(f"desc = {info['desc']}")
            lines.append("")

        # imports
        lines.append("import genice3.unitcell")
        lines.append("import numpy as np")
        lines.append("from genice2.cell import cellvectors")
        lines.append("")

        # クラス定義
        lines.append("")
        lines.append(f"class UnitCell(genice3.unitcell.UnitCell):")
        lines.append(f'    """')
        lines.append(f"    {info['module_name']}単位胞を定義するクラス。")
        lines.append("")
        lines.append("    NOTE: This unitcell is not yet implemented.")
        lines.append("    Please contact the maintainer or implement it manually.")
        lines.append(f'    """')
        lines.append("")
        lines.append("    def __init__(self, **kwargs):")
        lines.append("        raise NotImplementedError(")
        reason_str = ", ".join(reason) if reason else "自動変換できないため"
        lines.append(
            '            f"{self.__class__.__name__} is not yet implemented. "'
        )
        lines.append('            "This unitcell requires manual implementation. "')
        lines.append(
            '            "Please contact the maintainer or implement it manually. "'
        )
        lines.append(f'            f"Reason: {reason_str}"')
        lines.append("        )")

        return "\n".join(lines)

    # 通常の変換処理（以下は既存のコード）
    # docstringをそのまま移植
    if info.get("docstring"):
        lines.append(info["docstring"])
        lines.append("")

    # descをそのまま移植
    if info.get("desc"):
        lines.append(f"desc = {info['desc']}")
        lines.append("")

    # imports
    lines.append("import genice3.unitcell")
    lines.append("import numpy as np")
    if (
        info.get("pairs_str")
        or info.get("pairs")
        or info.get("fixed_str")
        or info.get("fixed")
    ):
        lines.append("import networkx as nx")
    lines.append("from genice2.cell import cellvectors")
    lines.append("")

    # クラス定義
    lines.append("")
    lines.append(f"class UnitCell(genice3.unitcell.UnitCell):")
    lines.append(f'    """')
    lines.append(f"    {info['module_name']}単位胞を定義するクラス。")
    lines.append(f'    """')
    lines.append("")
    lines.append("    def __init__(self, **kwargs):")

    # pairsの処理
    if info.get("pairs_str"):
        pairs = parse_pairs_string(info["pairs_str"])
        if pairs:
            lines.append('        pairs_str = """')
            for i, j in pairs:
                lines.append(f"        {i} {j}")
            lines.append('        """.split(')
            lines.append('            "\\n"')
            lines.append("        )")
            lines.append("        pairs = []")
            lines.append("        for line in pairs_str:")
            lines.append("            cols = line.split()")
            lines.append("            if len(cols) == 2:")
            lines.append("                pairs.append((int(cols[0]), int(cols[1])))")
            lines.append("")
            lines.append("        graph = nx.Graph(pairs)")
        else:
            lines.append("        graph = None  # pairsが空の場合は自動生成")
    elif info.get("pairs") == "fixed":
        # self.pairs = self.fixed の場合
        if info.get("fixed_str"):
            pairs = parse_pairs_string(info["fixed_str"])
            if pairs:
                lines.append('        pairs_str = """')
                for i, j in pairs:
                    lines.append(f"        {i} {j}")
                lines.append('        """.split(')
                lines.append('            "\\n"')
                lines.append("        )")
                lines.append("        pairs = []")
                lines.append("        for line in pairs_str:")
                lines.append("            cols = line.split()")
                lines.append("            if len(cols) == 2:")
                lines.append(
                    "                pairs.append((int(cols[0]), int(cols[1])))"
                )
                lines.append("")
                lines.append("        graph = nx.Graph(pairs)")
            else:
                lines.append("        graph = None")
        else:
            lines.append("        graph = None  # fixedから生成")
    elif info.get("has_cif"):
        lines.append("        # CIFから生成される場合はgraph=Noneで自動生成")
        lines.append("        graph = None")
    else:
        lines.append("        graph = None  # pairsがない場合は自動生成")

    # watersの処理
    if info.get("waters_str"):
        lines.append("")
        lines.append("        waters = np.fromstring(")
        lines.append('            """')
        # waters文字列を整形
        waters_lines = info["waters_str"].strip().split("\n")
        for line in waters_lines:
            if line.strip():
                lines.append(f"        {line.strip()}")
        lines.append('        """,')
        lines.append('            sep=" ",')
        lines.append("        ).reshape(-1, 3)")
    elif info.get("waters_array"):
        # np.array([...])形式の場合
        lines.append("")
        if isinstance(info["waters_array"], str):
            # ast.unparseで生成された文字列の場合
            # リストを整形して読みやすくする
            # 簡易的な整形: 各行が3要素のリストになるように改行を入れる
            array_str = info["waters_array"]
            # リストの各要素を抽出して整形
            # ただし、ast.unparseの結果をそのまま使う方が安全
            lines.append(f"        waters = np.array({array_str})")
        else:
            # ASTノードの場合（Python < 3.9）
            # 手動でリストを文字列に変換する必要がある
            # 簡易的な実装: 元のコードから抽出
            lines.append("        # np.array([...])形式のwaters")
            lines.append("        # 元のコードから手動で変換してください")
            lines.append("        waters = np.array([])  # TODO: 実装してください")

    # coord
    lines.append("")
    lines.append(f"        coord = \"{info['coord']}\"")

    # bondlen
    lines.append("")
    if isinstance(info["bondlen"], int):
        lines.append(f"        bondlen = {info['bondlen']}")
    else:
        lines.append(f"        bondlen = {info['bondlen']}")

    # density
    if info.get("density") is not None:
        lines.append("")
        lines.append(f"        # density = {info['density']}")

    # cell
    lines.append("")
    if info.get("cell"):
        # cellvectors()呼び出しがある場合
        lines.append(f"        cell = {info['cell']}")
    elif info.get("cell_str"):
        # 文字列からcellvectorsを生成（簡単な場合のみ）
        cell_vals = info["cell_str"].strip().split()
        if len(cell_vals) >= 3:
            a = cell_vals[0]
            b = cell_vals[1] if len(cell_vals) > 1 else a
            c = cell_vals[2] if len(cell_vals) > 2 else a
            lines.append(f"        cell = cellvectors(a={a}, b={b}, c={c})")
        else:
            lines.append(f"        # cell: {info['cell_str']}")
            lines.append("        # 手動でcellvectors()を設定してください")
    else:
        lines.append("        # cell: 手動で設定してください")

    # fixedの処理
    fixed_param = ""
    if info.get("fixed_str"):
        pairs = parse_pairs_string(info["fixed_str"])
        if pairs:
            lines.append("")
            lines.append("        fixed_pairs = [")
            for i, j in pairs:
                lines.append(f"            ({i}, {j}),")
            lines.append("        ]")
            lines.append("        fixed = nx.DiGraph(fixed_pairs)")
            fixed_param = ", fixed=fixed"
    elif info.get("fixed") == "pairs":
        lines.append("")
        lines.append("        # fixed = pairs の場合")
        lines.append("        fixed = nx.DiGraph(pairs) if pairs else nx.DiGraph()")
        fixed_param = ", fixed=fixed"

    # super().__init__呼び出し
    lines.append("")
    lines.append("        super().__init__(")
    lines.append("            cell=cell,")
    lines.append("            waters=waters,")
    # graphはNoneでない場合のみ渡す（Noneの場合は自動生成される）
    if info.get("pairs_str") or (
        info.get("pairs") == "fixed" and info.get("fixed_str")
    ):
        lines.append("            graph=graph,")
    lines.append("            coord=coord,")
    lines.append("            bondlen=bondlen,")
    if info.get("density") is not None:
        lines.append("            # density=density,")
    if fixed_param:
        lines.append(
            f"            {fixed_param[2:]},"
        )  # ', fixed=fixed' から ', ' を除去
    lines.append("            **kwargs,")
    lines.append("        )")

    return "\n".join(lines)


def convert_file(input_path: Path, output_dir: Optional[Path] = None) -> bool:
    """単一ファイルを変換

    自動変換できない場合は何も出力せずFalseを返す
    symlinkの場合はそのままコピーする
    """
    if not input_path.exists():
        return False

    # symlinkの場合はそのままコピー
    if input_path.is_symlink():
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            output_path = output_dir / f"{input_path.stem}.py"
        else:
            output_path = (
                input_path.parent.parent.parent
                / "genice3"
                / "unitcell"
                / f"{input_path.stem}.py"
            )
            output_path.parent.mkdir(parents=True, exist_ok=True)

        # symlinkのターゲットを取得（元のファイル名）
        target = input_path.readlink()

        # 相対パスに正規化
        if target.is_absolute():
            # 絶対パスの場合は、genice2/lattices/からの相対パスに変換
            lattices_dir = input_path.parent
            try:
                target = target.relative_to(lattices_dir)
            except ValueError:
                # 相対パスに変換できない場合はファイル名のみを使用
                target = Path(target.name)

        # 新しいディレクトリ内のファイルへのリンクにする
        # ターゲットファイル名のみを使用（同じディレクトリ内のファイルを指す）
        target_name = (
            target.name if hasattr(target, "name") else str(target).split("/")[-1]
        )
        new_target = Path(target_name)

        # 新しいsymlinkを作成
        if output_path.exists():
            output_path.unlink()
        output_path.symlink_to(new_target)
        return True

    info = analyze_lattice_file(input_path)

    if not info:
        return False

    code = generate_unitcell_code(info, generate_skeleton=True)

    # コードが生成できない場合は失敗
    if code is None:
        return False

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"{input_path.stem}.py"
    else:
        output_path = (
            input_path.parent.parent.parent
            / "genice3"
            / "unitcell"
            / f"{input_path.stem}.py"
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(code)

    return True


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    if sys.argv[1] == "--all":
        # すべてのlatticeファイルを変換
        lattices_dir = Path(__file__).parent.parent
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None

        lattice_files = sorted(lattices_dir.glob("*.py"))
        # __init__.pyや_toolディレクトリを除外
        lattice_files = [
            f
            for f in lattice_files
            if f.name != "__init__.py" and not f.name.startswith("__")
        ]

        success_count = 0
        for lattice_file in lattice_files:
            if convert_file(lattice_file, output_dir):
                success_count += 1

        print(f"\nConverted {success_count}/{len(lattice_files)} files")
        # すべて変換できた場合は成功、そうでない場合は失敗
        if success_count < len(lattice_files):
            sys.exit(1)
    else:
        # 単一ファイルを変換
        input_path = Path(sys.argv[1])
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None
        if not convert_file(input_path, output_dir):
            sys.exit(1)


if __name__ == "__main__":
    main()
