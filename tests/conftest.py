"""
Pytest configuration for all tests
プロジェクトのルートディレクトリをPythonパスに追加
"""

import sys
from pathlib import Path
import os

# プロジェクトのルートディレクトリをPythonパスに追加
# tests/conftest.py -> tests -> . と遡る
project_root = Path(__file__).resolve().parent.parent
project_root_str = str(project_root)

# 絶対パスで追加（確実にするため）
if project_root_str not in sys.path:
    sys.path.insert(0, project_root_str)

# 環境変数PYTHONPATHにも設定（pytestが読み込む前に）
if "PYTHONPATH" not in os.environ:
    os.environ["PYTHONPATH"] = project_root_str
elif project_root_str not in os.environ["PYTHONPATH"].split(os.pathsep):
    os.environ["PYTHONPATH"] = project_root_str + os.pathsep + os.environ["PYTHONPATH"]
