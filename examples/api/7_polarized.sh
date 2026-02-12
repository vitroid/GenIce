#!/bin/bash
# 分極した氷 (1) コンストラクタで target_pol を指定する場合のCLI（括弧付き）
# 対応するPython: examples/api/7_polarized_1.py

python3 -m genice3.cli.genice 1h \
  --rep 2 2 2 \
  --exporter _pol \
  --seed 114 \
  --depol_loop 1000 \
  --target_polarization 4 0 0
