"""
Command line: ./genice2.x 4 -f reshape[1,1,1,1,-1,0,1,1,-2]
Reshaping the unit cell.
  i:[1 1 1]
  j:[ 1 -1  0]
  k:[ 1  1 -2]

"""

from genice2.cell import cellvectors
import numpy as np
import genice2.lattices
desc = {"ref": {"IV": 'Avogadro'},
        "usage": "No options available.",
        "brief": "Orthogonalized ice IV."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.bondlen = 3.0000000000000004
        self.coord = 'relative'
        self.density = 1.3072141048893433
        self.waters = """
            0.2927    0.0000    1.0000
            0.4594    0.0000    0.8333
            0.4594    0.7500    0.0833
            0.6261    0.7500    0.9167
            0.4594    0.2500    0.0833
            0.6261    0.2500    0.9167
            0.6261    0.0000    0.1667
            0.7927    0.0000    1.0000
            0.4573    1.0000    1.0000
            0.6239    0.0000    0.8333
            0.6239    0.7500    0.0833
            0.7906    0.7500    0.9167
            0.6239    0.2500    0.0833
            0.7906    0.2500    0.9167
            0.7906    0.0000    0.1667
            0.9572    0.0000    1.0000
            0.2073    0.0000    0.0000
            0.3739    0.0000    0.8333
            0.3739    0.7500    0.0833
            0.5406    0.7500    0.9167
            0.3739    0.2500    0.0833
            0.5406    0.2500    0.9167
            0.5406    0.0000    0.1667
            0.7073    0.0000    1.0000
            0.0427    0.0000    0.0000
            0.2094    1.0000    0.8333
            0.2094    0.7500    0.0833
            0.3761    0.7500    0.9167
            0.2094    0.2500    0.0833
            0.3761    0.2500    0.9167
            0.3761    1.0000    0.1667
            0.5427    0.0000    1.0000
            0.2550    0.9678    0.9074
            0.4216    0.9678    0.7407
            0.4216    0.7178    0.9907
            0.5883    0.7178    0.8241
            0.4216    0.2178    0.9907
            0.5883    0.2178    0.8241
            0.5883    0.9678    0.0741
            0.7550    0.9678    0.9074
            0.1617    0.9678    0.9259
            0.3284    0.9678    0.7593
            0.3284    0.7178    0.0093
            0.4950    0.7178    0.8426
            0.3284    0.2178    0.0093
            0.4950    0.2178    0.8426
            0.4950    0.9678    0.0926
            0.6617    0.9678    0.9259
            0.2450    0.1228    0.9376
            0.4117    0.1228    0.7709
            0.4117    0.8728    0.0209
            0.5784    0.8728    0.8543
            0.4117    0.3728    0.0209
            0.5783    0.3728    0.8543
            0.5783    0.1228    0.1043
            0.7450    0.1228    0.9376
            0.3383    0.0950    0.9469
            0.5050    0.0950    0.7802
            0.5050    0.8450    0.0302
            0.6716    0.8450    0.8636
            0.5050    0.3450    0.0302
            0.6716    0.3450    0.8636
            0.6716    0.0950    0.1136
            0.8383    0.0950    0.9469
            0.2550    0.1550    0.0302
            0.4216    0.1550    0.8636
            0.4216    0.9050    0.1136
            0.5883    0.9050    0.9469
            0.4216    0.4050    0.1136
            0.5883    0.4050    0.9469
            0.5883    0.1550    0.1969
            0.7550    0.1550    0.0302
            0.1617    0.1272    0.0209
            0.3284    0.1272    0.8543
            0.3283    0.8772    0.1043
            0.4950    0.8772    0.9376
            0.3283    0.3772    0.1043
            0.4950    0.3772    0.9376
            0.4950    0.1272    0.1876
            0.6617    0.1272    0.0209
            0.2450    0.0322    0.0926
            0.4117    0.0322    0.9259
            0.4117    0.7822    0.1759
            0.5784    0.7822    0.0093
            0.4117    0.2822    0.1759
            0.5784    0.2822    0.0093
            0.5783    0.0322    0.2593
            0.7450    0.0322    0.0926
            0.3383    0.0322    0.0741
            0.5050    0.0322    0.9074
            0.5050    0.7822    0.1574
            0.6716    0.7822    0.9907
            0.5050    0.2822    0.1574
            0.6716    0.2822    0.9907
            0.6716    0.0322    0.2407
            0.8383    0.0322    0.0741
            0.2550    0.8772    0.0624
            0.4216    0.8772    0.8957
            0.4216    0.6272    0.1457
            0.5883    0.6272    0.9791
            0.4216    0.1272    0.1457
            0.5883    0.1272    0.9791
            0.5883    0.8772    0.2291
            0.7550    0.8772    0.0624
            0.1617    0.9050    0.0531
            0.3284    0.9050    0.8864
            0.3283    0.6550    0.1364
            0.4950    0.6550    0.9698
            0.3284    0.1550    0.1364
            0.4950    0.1550    0.9698
            0.4950    0.9050    0.2198
            0.6617    0.9050    0.0531
            0.2450    0.8450    0.9698
            0.4117    0.8450    0.8031
            0.4117    0.5950    0.0531
            0.5783    0.5950    0.8864
            0.4117    0.0950    0.0531
            0.5783    0.0950    0.8864
            0.5784    0.8450    0.1364
            0.7450    0.8450    0.9698
            0.3383    0.8728    0.9791
            0.5050    0.8728    0.8124
            0.5050    0.6228    0.0624
            0.6716    0.6228    0.8957
            0.5050    0.1228    0.0624
            0.6716    0.1228    0.8957
            0.6716    0.8728    0.1457
            0.8383    0.8728    0.9791
            0.6261    0.0000    0.6667
            0.7928    0.0000    0.5000
            0.7927    0.7500    0.7500
            0.9594    0.7500    0.5833
            0.7927    0.2500    0.7500
            0.9594    0.2500    0.5833
            0.9594    1.0000    0.8333
            0.1261    0.0000    0.6667
            0.7906    1.0000    0.6667
            0.9573    1.0000    0.5000
            0.9573    0.7500    0.7500
            0.1239    0.7500    0.5833
            0.9573    0.2500    0.7500
            0.1239    0.2500    0.5833
            0.1239    1.0000    0.8333
            0.2906    1.0000    0.6667
            0.5406    0.0000    0.6667
            0.7073    0.0000    0.5000
            0.7073    0.7500    0.7500
            0.8739    0.7500    0.5833
            0.7073    0.2500    0.7500
            0.8739    0.2500    0.5833
            0.8739    0.0000    0.8333
            0.0406    0.0000    0.6667
            0.3761    1.0000    0.6667
            0.5427    0.0000    0.5000
            0.5428    0.7500    0.7500
            0.7094    0.7500    0.5833
            0.5427    0.2500    0.7500
            0.7094    0.2500    0.5833
            0.7094    1.0000    0.8333
            0.8761    0.0000    0.6667
            0.5883    0.9678    0.5741
            0.7550    0.9678    0.4074
            0.7550    0.7178    0.6574
            0.9217    0.7178    0.4907
            0.7550    0.2178    0.6574
            0.9217    0.2178    0.4907
            0.9216    0.9678    0.7407
            0.0883    0.9678    0.5741
            0.4950    0.9678    0.5926
            0.6617    0.9678    0.4259
            0.6617    0.7178    0.6759
            0.8284    0.7178    0.5093
            0.6617    0.2178    0.6759
            0.8284    0.2178    0.5093
            0.8284    0.9678    0.7593
            0.9950    0.9678    0.5926
            0.5784    0.1228    0.6043
            0.7450    0.1228    0.4376
            0.7450    0.8728    0.6876
            0.9117    0.8728    0.5209
            0.7450    0.3728    0.6876
            0.9117    0.3728    0.5209
            0.9117    0.1228    0.7709
            0.0784    0.1228    0.6043
            0.6717    0.0950    0.6136
            0.8383    0.0950    0.4469
            0.8383    0.8450    0.6969
            0.0050    0.8450    0.5302
            0.8383    0.3450    0.6969
            0.0050    0.3450    0.5302
            0.0050    0.0950    0.7802
            0.1716    0.0950    0.6136
            0.5883    0.1550    0.6969
            0.7550    0.1550    0.5302
            0.7550    0.9050    0.7802
            0.9216    0.9050    0.6136
            0.7550    0.4050    0.7802
            0.9216    0.4050    0.6136
            0.9216    0.1550    0.8636
            0.0883    0.1550    0.6969
            0.4950    0.1272    0.6876
            0.6617    0.1272    0.5209
            0.6617    0.8772    0.7709
            0.8284    0.8772    0.6043
            0.6617    0.3772    0.7709
            0.8284    0.3772    0.6043
            0.8284    0.1272    0.8543
            0.9950    0.1272    0.6876
            0.5784    0.0322    0.7593
            0.7450    0.0322    0.5926
            0.7450    0.7822    0.8426
            0.9117    0.7822    0.6759
            0.7450    0.2822    0.8426
            0.9117    0.2822    0.6759
            0.9117    0.0322    0.9259
            0.0784    0.0322    0.7593
            0.6716    0.0322    0.7407
            0.8383    0.0322    0.5741
            0.8383    0.7822    0.8241
            0.0050    0.7822    0.6574
            0.8383    0.2822    0.8241
            0.0050    0.2822    0.6574
            0.0050    0.0322    0.9074
            0.1716    0.0322    0.7407
            0.5883    0.8772    0.7291
            0.7550    0.8772    0.5624
            0.7550    0.6272    0.8124
            0.9216    0.6272    0.6457
            0.7550    0.1272    0.8124
            0.9216    0.1272    0.6457
            0.9216    0.8772    0.8957
            0.0883    0.8772    0.7291
            0.4950    0.9050    0.7198
            0.6617    0.9050    0.5531
            0.6617    0.6550    0.8031
            0.8284    0.6550    0.6364
            0.6617    0.1550    0.8031
            0.8284    0.1550    0.6364
            0.8284    0.9050    0.8864
            0.9950    0.9050    0.7198
            0.5784    0.8450    0.6364
            0.7450    0.8450    0.4698
            0.7450    0.5950    0.7198
            0.9117    0.5950    0.5531
            0.7450    0.0950    0.7198
            0.9117    0.0950    0.5531
            0.9117    0.8450    0.8031
            0.0783    0.8450    0.6364
            0.6716    0.8728    0.6457
            0.8383    0.8728    0.4791
            0.8383    0.6228    0.7291
            0.0050    0.6228    0.5624
            0.8383    0.1228    0.7291
            0.0050    0.1228    0.5624
            0.0050    0.8728    0.8124
            0.1716    0.8728    0.6457
            0.9594    0.0000    0.3333
            0.1261    0.0000    0.1667
            0.1261    0.7500    0.4167
            0.2927    0.7500    0.2500
            0.1261    0.2500    0.4167
            0.2927    0.2500    0.2500
            0.2927    1.0000    0.5000
            0.4594    0.0000    0.3333
            0.1239    1.0000    0.3333
            0.2906    0.0000    0.1667
            0.2906    0.7500    0.4167
            0.4573    0.7500    0.2500
            0.2906    0.2500    0.4167
            0.4573    0.2500    0.2500
            0.4573    1.0000    0.5000
            0.6239    1.0000    0.3333
            0.8739    0.0000    0.3333
            0.0406    0.0000    0.1667
            0.0406    0.7500    0.4167
            0.2073    0.7500    0.2500
            0.0406    0.2500    0.4167
            0.2073    0.2500    0.2500
            0.2073    0.0000    0.5000
            0.3739    0.0000    0.3333
            0.7094    0.0000    0.3333
            0.8761    0.0000    0.1667
            0.8761    0.7500    0.4167
            0.0427    0.7500    0.2500
            0.8761    0.2500    0.4167
            0.0427    0.2500    0.2500
            0.0427    1.0000    0.5000
            0.2094    0.0000    0.3333
            0.9216    0.9678    0.2407
            0.0883    0.9678    0.0741
            0.0883    0.7178    0.3241
            0.2550    0.7178    0.1574
            0.0883    0.2178    0.3241
            0.2550    0.2178    0.1574
            0.2550    0.9678    0.4074
            0.4217    0.9678    0.2407
            0.8284    0.9678    0.2593
            0.9950    0.9678    0.0926
            0.9950    0.7178    0.3426
            0.1617    0.7178    0.1759
            0.9950    0.2178    0.3426
            0.1617    0.2178    0.1759
            0.1617    0.9678    0.4259
            0.3284    0.9678    0.2593
            0.9117    0.1228    0.2709
            0.0784    0.1228    0.1043
            0.0784    0.8728    0.3543
            0.2450    0.8728    0.1876
            0.0784    0.3728    0.3543
            0.2450    0.3728    0.1876
            0.2450    0.1228    0.4376
            0.4117    0.1228    0.2709
            0.0050    0.0950    0.2802
            0.1717    0.0950    0.1136
            0.1716    0.8450    0.3636
            0.3383    0.8450    0.1969
            0.1716    0.3450    0.3636
            0.3383    0.3450    0.1969
            0.3383    0.0950    0.4469
            0.5050    0.0950    0.2802
            0.9216    0.1550    0.3636
            0.0883    0.1550    0.1969
            0.0883    0.9050    0.4469
            0.2550    0.9050    0.2802
            0.0883    0.4050    0.4469
            0.2550    0.4050    0.2802
            0.2550    0.1550    0.5302
            0.4216    0.1550    0.3636
            0.8284    0.1272    0.3543
            0.9950    0.1272    0.1876
            0.9950    0.8772    0.4376
            0.1617    0.8772    0.2709
            0.9950    0.3772    0.4376
            0.1617    0.3772    0.2709
            0.1617    0.1272    0.5209
            0.3284    0.1272    0.3543
            0.9117    0.0322    0.4259
            0.0784    0.0322    0.2593
            0.0784    0.7822    0.5093
            0.2450    0.7822    0.3426
            0.0784    0.2822    0.5093
            0.2450    0.2822    0.3426
            0.2450    0.0322    0.5926
            0.4117    0.0322    0.4259
            0.0050    0.0322    0.4074
            0.1716    0.0322    0.2407
            0.1716    0.7822    0.4907
            0.3383    0.7822    0.3241
            0.1716    0.2822    0.4907
            0.3383    0.2822    0.3241
            0.3383    0.0322    0.5741
            0.5050    0.0322    0.4074
            0.9217    0.8772    0.3957
            0.0883    0.8772    0.2291
            0.0883    0.6272    0.4791
            0.2550    0.6272    0.3124
            0.0883    0.1272    0.4791
            0.2550    0.1272    0.3124
            0.2550    0.8772    0.5624
            0.4216    0.8772    0.3957
            0.8284    0.9050    0.3864
            0.9950    0.9050    0.2198
            0.9950    0.6550    0.4698
            0.1617    0.6550    0.3031
            0.9950    0.1550    0.4698
            0.1617    0.1550    0.3031
            0.1617    0.9050    0.5531
            0.3284    0.9050    0.3864
            0.9117    0.8450    0.3031
            0.0784    0.8450    0.1364
            0.0784    0.5950    0.3864
            0.2450    0.5950    0.2198
            0.0784    0.0950    0.3864
            0.2450    0.0950    0.2198
            0.2450    0.8450    0.4698
            0.4117    0.8450    0.3031
            0.0050    0.8728    0.3124
            0.1716    0.8728    0.1457
            0.1716    0.6228    0.3957
            0.3383    0.6228    0.2291
            0.1716    0.1228    0.3957
            0.3383    0.1228    0.2291
            0.3383    0.8728    0.4791
            0.5050    0.8728    0.3124
            0.6261    0.5000    0.1667
            0.7927    0.5000    1.0000
            0.7927    0.2500    0.2500
            0.9594    0.2500    0.0833
            0.7927    0.7500    0.2500
            0.9594    0.7500    0.0833
            0.9594    0.5000    0.3333
            0.1261    0.5000    0.1667
            0.7906    0.5000    0.1667
            0.9573    0.5000    1.0000
            0.9573    0.2500    0.2500
            0.1239    0.2500    0.0833
            0.9573    0.7500    0.2500
            0.1239    0.7500    0.0833
            0.1239    0.5000    0.3333
            0.2906    0.5000    0.1667
            0.5406    0.5000    0.1667
            0.7073    0.5000    0.0000
            0.7073    0.2500    0.2500
            0.8739    0.2500    0.0833
            0.7073    0.7500    0.2500
            0.8739    0.7500    0.0833
            0.8739    0.5000    0.3333
            0.0406    0.5000    0.1667
            0.3761    0.5000    0.1667
            0.5427    0.5000    1.0000
            0.5427    0.2500    0.2500
            0.7094    0.2500    0.0833
            0.5427    0.7500    0.2500
            0.7094    0.7500    0.0833
            0.7094    0.5000    0.3333
            0.8761    0.5000    0.1667
            0.5883    0.4678    0.0741
            0.7550    0.4678    0.9074
            0.7550    0.2178    0.1574
            0.9217    0.2178    0.9907
            0.7550    0.7178    0.1574
            0.9217    0.7178    0.9907
            0.9216    0.4678    0.2407
            0.0883    0.4678    0.0741
            0.4950    0.4678    0.0926
            0.6617    0.4678    0.9259
            0.6617    0.2178    0.1759
            0.8284    0.2178    0.0093
            0.6617    0.7178    0.1759
            0.8284    0.7178    0.0093
            0.8284    0.4678    0.2593
            0.9950    0.4678    0.0926
            0.5784    0.6228    0.1043
            0.7450    0.6228    0.9376
            0.7450    0.3728    0.1876
            0.9117    0.3728    0.0209
            0.7450    0.8728    0.1876
            0.9117    0.8728    0.0209
            0.9117    0.6228    0.2709
            0.0783    0.6228    0.1043
            0.6716    0.5950    0.1136
            0.8383    0.5950    0.9469
            0.8383    0.3450    0.1969
            0.0050    0.3450    0.0302
            0.8383    0.8450    0.1969
            0.0050    0.8450    0.0302
            0.0050    0.5950    0.2802
            0.1716    0.5950    0.1136
            0.5883    0.6550    0.1969
            0.7550    0.6550    0.0302
            0.7550    0.4050    0.2802
            0.9216    0.4050    0.1136
            0.7550    0.9050    0.2802
            0.9216    0.9050    0.1136
            0.9216    0.6550    0.3636
            0.0883    0.6550    0.1969
            0.4950    0.6272    0.1876
            0.6617    0.6272    0.0209
            0.6617    0.3772    0.2709
            0.8284    0.3772    0.1043
            0.6617    0.8772    0.2709
            0.8284    0.8772    0.1043
            0.8283    0.6272    0.3543
            0.9950    0.6272    0.1876
            0.5784    0.5322    0.2593
            0.7450    0.5322    0.0926
            0.7450    0.2822    0.3426
            0.9117    0.2822    0.1759
            0.7450    0.7822    0.3426
            0.9117    0.7822    0.1759
            0.9117    0.5322    0.4259
            0.0784    0.5322    0.2593
            0.6717    0.5322    0.2407
            0.8383    0.5322    0.0741
            0.8383    0.2822    0.3241
            0.0050    0.2822    0.1574
            0.8383    0.7822    0.3241
            0.0050    0.7822    0.1574
            0.0050    0.5322    0.4074
            0.1716    0.5322    0.2407
            0.5883    0.3772    0.2291
            0.7550    0.3772    0.0624
            0.7550    0.1272    0.3124
            0.9216    0.1272    0.1457
            0.7550    0.6272    0.3124
            0.9216    0.6272    0.1457
            0.9216    0.3772    0.3957
            0.0883    0.3772    0.2291
            0.4950    0.4050    0.2198
            0.6617    0.4050    0.0531
            0.6617    0.1550    0.3031
            0.8284    0.1550    0.1364
            0.6617    0.6550    0.3031
            0.8284    0.6550    0.1364
            0.8284    0.4050    0.3864
            0.9950    0.4050    0.2198
            0.5784    0.3450    0.1364
            0.7450    0.3450    0.9698
            0.7450    0.0950    0.2198
            0.9117    0.0950    0.0531
            0.7450    0.5950    0.2198
            0.9117    0.5950    0.0531
            0.9117    0.3450    0.3031
            0.0783    0.3450    0.1364
            0.6717    0.3728    0.1457
            0.8383    0.3728    0.9791
            0.8383    0.1228    0.2291
            0.0050    0.1228    0.0624
            0.8383    0.6228    0.2291
            0.0050    0.6228    0.0624
            0.0050    0.3728    0.3124
            0.1716    0.3728    0.1457
            0.9594    0.5000    0.8333
            0.1261    0.5000    0.6667
            0.1261    0.2500    0.9167
            0.2927    0.2500    0.7500
            0.1261    0.7500    0.9167
            0.2927    0.7500    0.7500
            0.2927    0.5000    1.0000
            0.4594    0.5000    0.8333
            0.1239    0.5000    0.8333
            0.2906    0.5000    0.6667
            0.2906    0.2500    0.9167
            0.4573    0.2500    0.7500
            0.2906    0.7500    0.9167
            0.4573    0.7500    0.7500
            0.4573    0.5000    1.0000
            0.6239    0.5000    0.8333
            0.8739    0.5000    0.8333
            0.0406    0.5000    0.6667
            0.0406    0.2500    0.9167
            0.2073    0.2500    0.7500
            0.0406    0.7500    0.9167
            0.2073    0.7500    0.7500
            0.2073    0.5000    1.0000
            0.3739    0.5000    0.8333
            0.7094    0.5000    0.8333
            0.8761    0.5000    0.6667
            0.8761    0.2500    0.9167
            0.0427    0.2500    0.7500
            0.8761    0.7500    0.9167
            0.0427    0.7500    0.7500
            0.0427    0.5000    1.0000
            0.2094    0.5000    0.8333
            0.9216    0.4678    0.7407
            0.0883    0.4678    0.5741
            0.0883    0.2178    0.8241
            0.2550    0.2178    0.6574
            0.0883    0.7178    0.8241
            0.2550    0.7178    0.6574
            0.2550    0.4678    0.9074
            0.4216    0.4678    0.7407
            0.8284    0.4678    0.7593
            0.9950    0.4678    0.5926
            0.9950    0.2178    0.8426
            0.1617    0.2178    0.6759
            0.9950    0.7178    0.8426
            0.1617    0.7178    0.6759
            0.1617    0.4678    0.9259
            0.3284    0.4678    0.7593
            0.9117    0.6228    0.7709
            0.0784    0.6228    0.6043
            0.0784    0.3728    0.8543
            0.2450    0.3728    0.6876
            0.0784    0.8728    0.8543
            0.2450    0.8728    0.6876
            0.2450    0.6228    0.9376
            0.4117    0.6228    0.7709
            0.0050    0.5950    0.7802
            0.1716    0.5950    0.6136
            0.1716    0.3450    0.8636
            0.3383    0.3450    0.6969
            0.1716    0.8450    0.8636
            0.3383    0.8450    0.6969
            0.3383    0.5950    0.9469
            0.5050    0.5950    0.7802
            0.9216    0.6550    0.8636
            0.0883    0.6550    0.6969
            0.0883    0.4050    0.9469
            0.2550    0.4050    0.7802
            0.0883    0.9050    0.9469
            0.2550    0.9050    0.7802
            0.2550    0.6550    0.0302
            0.4216    0.6550    0.8636
            0.8284    0.6272    0.8543
            0.9950    0.6272    0.6876
            0.9950    0.3772    0.9376
            0.1617    0.3772    0.7709
            0.9950    0.8772    0.9376
            0.1617    0.8772    0.7709
            0.1617    0.6272    0.0209
            0.3284    0.6272    0.8543
            0.9117    0.5322    0.9259
            0.0784    0.5322    0.7593
            0.0784    0.2822    0.0093
            0.2450    0.2822    0.8426
            0.0784    0.7822    0.0093
            0.2450    0.7822    0.8426
            0.2450    0.5322    0.0926
            0.4117    0.5322    0.9259
            0.0050    0.5322    0.9074
            0.1716    0.5322    0.7407
            0.1716    0.2822    0.9907
            0.3383    0.2822    0.8241
            0.1716    0.7822    0.9907
            0.3383    0.7822    0.8241
            0.3383    0.5322    0.0741
            0.5050    0.5322    0.9074
            0.9216    0.3772    0.8957
            0.0883    0.3772    0.7291
            0.0883    0.1272    0.9791
            0.2550    0.1272    0.8124
            0.0883    0.6272    0.9791
            0.2550    0.6272    0.8124
            0.2550    0.3772    0.0624
            0.4216    0.3772    0.8957
            0.8284    0.4050    0.8864
            0.9950    0.4050    0.7198
            0.9950    0.1550    0.9698
            0.1617    0.1550    0.8031
            0.9950    0.6550    0.9698
            0.1617    0.6550    0.8031
            0.1617    0.4050    0.0531
            0.3284    0.4050    0.8864
            0.9117    0.3450    0.8031
            0.0784    0.3450    0.6364
            0.0784    0.0950    0.8864
            0.2450    0.0950    0.7198
            0.0784    0.5950    0.8864
            0.2450    0.5950    0.7198
            0.2450    0.3450    0.9698
            0.4117    0.3450    0.8031
            0.0050    0.3728    0.8124
            0.1716    0.3728    0.6457
            0.1716    0.1228    0.8957
            0.3383    0.1228    0.7291
            0.1716    0.6228    0.8957
            0.3383    0.6228    0.7291
            0.3383    0.3728    0.9791
            0.5050    0.3728    0.8124
            0.2927    0.5000    0.5000
            0.4594    0.5000    0.3333
            0.4594    0.2500    0.5833
            0.6261    0.2500    0.4167
            0.4594    0.7500    0.5833
            0.6261    0.7500    0.4167
            0.6261    0.5000    0.6667
            0.7927    0.5000    0.5000
            0.4573    0.5000    0.5000
            0.6239    0.5000    0.3333
            0.6239    0.2500    0.5833
            0.7906    0.2500    0.4167
            0.6239    0.7500    0.5833
            0.7906    0.7500    0.4167
            0.7906    0.5000    0.6667
            0.9573    0.5000    0.5000
            0.2073    0.5000    0.5000
            0.3739    0.5000    0.3333
            0.3739    0.2500    0.5833
            0.5406    0.2500    0.4167
            0.3739    0.7500    0.5833
            0.5406    0.7500    0.4167
            0.5406    0.5000    0.6667
            0.7073    0.5000    0.5000
            0.0427    0.5000    0.5000
            0.2094    0.5000    0.3333
            0.2094    0.2500    0.5833
            0.3761    0.2500    0.4167
            0.2094    0.7500    0.5833
            0.3761    0.7500    0.4167
            0.3761    0.5000    0.6667
            0.5427    0.5000    0.5000
            0.2550    0.4678    0.4074
            0.4217    0.4678    0.2407
            0.4216    0.2178    0.4907
            0.5883    0.2178    0.3241
            0.4216    0.7178    0.4907
            0.5883    0.7178    0.3241
            0.5883    0.4678    0.5741
            0.7550    0.4678    0.4074
            0.1617    0.4678    0.4259
            0.3284    0.4678    0.2593
            0.3284    0.2178    0.5093
            0.4950    0.2178    0.3426
            0.3284    0.7178    0.5093
            0.4950    0.7178    0.3426
            0.4950    0.4678    0.5926
            0.6617    0.4678    0.4259
            0.2450    0.6228    0.4376
            0.4117    0.6228    0.2709
            0.4117    0.3728    0.5209
            0.5784    0.3728    0.3543
            0.4117    0.8728    0.5209
            0.5784    0.8728    0.3543
            0.5783    0.6228    0.6043
            0.7450    0.6228    0.4376
            0.3383    0.5950    0.4469
            0.5050    0.5950    0.2802
            0.5050    0.3450    0.5302
            0.6716    0.3450    0.3636
            0.5050    0.8450    0.5302
            0.6716    0.8450    0.3636
            0.6716    0.5950    0.6136
            0.8383    0.5950    0.4469
            0.2550    0.6550    0.5302
            0.4217    0.6550    0.3636
            0.4216    0.4050    0.6136
            0.5883    0.4050    0.4469
            0.4216    0.9050    0.6136
            0.5883    0.9050    0.4469
            0.5883    0.6550    0.6969
            0.7550    0.6550    0.5302
            0.1617    0.6272    0.5209
            0.3284    0.6272    0.3543
            0.3284    0.3772    0.6043
            0.4950    0.3772    0.4376
            0.3284    0.8772    0.6043
            0.4950    0.8772    0.4376
            0.4950    0.6272    0.6876
            0.6617    0.6272    0.5209
            0.2450    0.5322    0.5926
            0.4117    0.5322    0.4259
            0.4117    0.2822    0.6759
            0.5784    0.2822    0.5093
            0.4117    0.7822    0.6759
            0.5783    0.7822    0.5093
            0.5783    0.5322    0.7593
            0.7450    0.5322    0.5926
            0.3383    0.5322    0.5741
            0.5050    0.5322    0.4074
            0.5050    0.2822    0.6574
            0.6716    0.2822    0.4907
            0.5050    0.7822    0.6574
            0.6716    0.7822    0.4907
            0.6716    0.5322    0.7407
            0.8383    0.5322    0.5741
            0.2550    0.3772    0.5624
            0.4216    0.3772    0.3957
            0.4216    0.1272    0.6457
            0.5883    0.1272    0.4791
            0.4216    0.6272    0.6457
            0.5883    0.6272    0.4791
            0.5883    0.3772    0.7291
            0.7550    0.3772    0.5624
            0.1617    0.4050    0.5531
            0.3284    0.4050    0.3864
            0.3284    0.1550    0.6364
            0.4950    0.1550    0.4698
            0.3284    0.6550    0.6364
            0.4950    0.6550    0.4698
            0.4950    0.4050    0.7198
            0.6617    0.4050    0.5531
            0.2450    0.3450    0.4698
            0.4117    0.3450    0.3031
            0.4117    0.0950    0.5531
            0.5784    0.0950    0.3864
            0.4117    0.5950    0.5531
            0.5783    0.5950    0.3864
            0.5784    0.3450    0.6364
            0.7450    0.3450    0.4698
            0.3383    0.3728    0.4791
            0.5050    0.3728    0.3124
            0.5050    0.1228    0.5624
            0.6716    0.1228    0.3957
            0.5050    0.6228    0.5624
            0.6716    0.6228    0.3957
            0.6716    0.3728    0.6457
            0.8383    0.3728    0.4791
        """

        self.waters = np.fromstring(self.waters, dtype='float', sep=" ")
        self.waters = self.waters.reshape(self.waters.shape[0] // 3, 3)
        self.waters -= np.floor(self.waters)

        halfwater = []
        for w in self.waters:
            if np.all(w < 0.5):
                halfwater.append(w)
        self.waters = np.array(halfwater) * 2

        self.cell = cellvectors(a=16.919963345,
                                b=8.65462208,
                                c=14.99024516)
