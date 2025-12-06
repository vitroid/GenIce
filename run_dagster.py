#!/usr/bin/env python3
"""
Dagster UIを起動するスクリプト。

使用方法:
    python run_dagster.py

または:
    dagster dev -m genice2.genice3_dagster
"""

import sys
from dagster import DagsterInstance
from dagster._cli import main as dagster_main

if __name__ == "__main__":
    # Dagsterのコマンドライン引数を設定
    sys.argv = ["dagster", "dev", "-m", "genice2.genice3_dagster"]
    dagster_main()
