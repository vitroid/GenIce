#!/bin/bash
# GenIce3の最大限のオプションを含む例
# 対応するPythonコード: examples/api/3_doped.py
# 対応するYAML設定ファイル: examples/api/3_doped.yaml

python3 -m genice3.cli.genice A15 --shift 0.1 0.1 0.1 --anion 15=Cl --cation 21=Na --density 0.8 \
  --rep 2 2 2 \
  --exporter gromacs --water_model 4site \
  --seed 42 \
  --depol_loop 2000 \
  --spot_anion 1=Cl \
  --spot_cation 5=Na
