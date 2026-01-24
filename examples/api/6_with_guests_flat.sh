#!/bin/bash
# GenIce3の最大限のオプションを含む例
# 対応するPythonコード: examples/api/3.py
# 対応するTOML設定ファイル: examples/api/3.toml

python3 -m genice3.cli.genice A15 --shift 0.1 0.1 0.1 --density 0.8 \
  --rep 2 2 2 \
  --exporter gromacs --guest A12=me --guest A14=et --spot_guest 0=4site --water_model 4site \
  --seed 42 \
  --depol_loop 2000
