# coding: utf-8
desc={"ref": {"II+IVa":         'Karttunen 2011',
              "T":              'Sikiric 2010'},
      "usage": "No options available.",
      "brief": "Hypothetical clathrate type T."
      }

import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    """
[T] A. J. Karttunen, T. F. Fässler, M. Linnolahti, T. A. Pakkanen, Inorg Chem, 2011, DOI:10.1021/ic102178d.
"""
    def __init__(self):
        self.cell="""
        2.75276384094 2.75276384094 2.75276384094
        """

        self.bondlen=0.28
        self.coord="relative"

        self.waters="""
        0.062500 0.937500 0.937500
        0.961373 0.898873 1.000000
        0.937500 0.937500 0.937500
        0.038627 0.898873 1.000000
        1.000000 0.961373 0.898873
        0.062500 0.062500 0.937500
        0.038627 0.101127 1.000000
        0.961373 0.101127 1.000000
        0.937500 0.937500 0.062500
        0.898873 1.000000 0.961373
        0.898873 1.000000 0.038627
        0.937500 0.062500 0.937500
        1.000000 0.038627 0.898873
        0.062500 0.062500 0.062500
        0.062500 0.937500 0.062500
        1.000000 0.961373 0.101127
        0.937500 0.062500 0.062500
        1.000000 0.038627 0.101127
        0.101127 1.000000 0.961373
        0.101127 1.000000 0.038627
        0.078125 0.759119 0.881506
        1.000000 0.727869 0.928381
        0.046875 0.805994 0.803381
        0.953125 0.805994 0.803381
        0.921875 0.759119 0.881506
        0.062500 0.619364 0.740011
        0.109375 0.697489 0.771261
        0.070748 0.767367 0.740881
        0.069877 0.618494 0.873258
        0.038627 0.571619 0.795133
        0.116752 0.696619 0.857633
        0.937500 0.619364 0.740011
        0.890625 0.697489 0.771261
        1.000000 0.643237 0.701383
        1.000000 0.736987 0.701383
        0.929252 0.767367 0.740881
        0.961373 0.571619 0.795133
        0.930123 0.618494 0.873258
        0.883248 0.696619 0.857633
        0.000000 0.641496 0.921004
        0.883248 0.812500 0.929252
        0.930123 0.820748 1.000000
        0.890625 0.890625 0.890625
        0.929252 0.883248 0.812500
        0.000000 0.930123 0.820748
        0.069877 0.820748 1.000000
        1.000000 0.766496 0.000000
        0.109375 0.890625 0.890625
        0.116752 0.812500 0.929252
        0.070748 0.883248 0.812500
        0.767367 0.740881 0.070748
        0.619364 0.740011 0.062500
        0.697489 0.771261 0.109375
        0.571619 0.795133 0.038627
        0.618494 0.873258 0.069877
        0.696619 0.857633 0.116752
        0.805994 0.803381 0.046875
        0.805994 0.803381 0.953125
        0.727869 0.928381 1.000000
        0.759119 0.881506 0.078125
        0.759119 0.881506 0.921875
        0.641496 0.921004 1.000000
        0.643237 0.701383 0.000000
        0.736987 0.701383 1.000000
        0.619364 0.740011 0.937500
        0.697489 0.771261 0.890625
        0.767367 0.740881 0.929252
        0.571619 0.795133 0.961373
        0.618494 0.873258 0.930123
        0.696619 0.857633 0.883248
        0.766496 0.000000 1.000000
        0.812500 0.929252 0.883248
        0.820748 1.000000 0.930123
        0.883248 0.812500 0.070748
        0.890625 0.890625 0.109375
        0.820748 1.000000 0.069877
        0.812500 0.929252 0.116752
        0.881506 0.078125 0.759119
        0.928381 1.000000 0.727869
        0.803381 0.046875 0.805994
        0.803381 0.953125 0.805994
        0.881506 0.921875 0.759119
        0.740881 0.070748 0.767367
        0.740011 0.062500 0.619364
        0.771261 0.109375 0.697489
        0.857633 0.116752 0.696619
        0.873258 0.069877 0.618494
        0.795133 0.038627 0.571619
        0.701383 0.000000 0.643237
        0.701383 1.000000 0.736987
        0.740881 0.929252 0.767367
        0.740011 0.937500 0.619364
        0.771261 0.890625 0.697489
        0.795133 0.961373 0.571619
        0.873258 0.930123 0.618494
        0.857633 0.883248 0.696619
        0.921004 1.000000 0.641496
        0.929252 0.116752 0.812500
        0.000000 0.069877 0.820748
        0.890625 0.109375 0.890625
        0.812500 0.070748 0.883248
        0.000000 1.000000 0.766496
        0.437500 0.437500 0.562500
        0.461373 0.398873 0.500000
        0.538627 0.398873 0.500000
        0.538627 0.601127 0.500000
        0.562500 0.562500 0.562500
        0.500000 0.538627 0.601127
        0.461373 0.601127 0.500000
        0.437500 0.562500 0.562500
        0.601127 0.500000 0.461373
        0.601127 0.500000 0.538627
        0.562500 0.562500 0.437500
        0.500000 0.461373 0.601127
        0.562500 0.437500 0.562500
        0.437500 0.437500 0.437500
        0.500000 0.538627 0.398873
        0.437500 0.562500 0.437500
        0.562500 0.437500 0.437500
        0.500000 0.461373 0.398873
        0.398873 0.500000 0.538627
        0.398873 0.500000 0.461373
        0.383248 0.803381 0.642367
        0.430123 0.881506 0.616325
        0.461373 0.928381 0.699658
        0.390625 0.802511 0.728739
        0.437500 0.880636 0.765198
        0.500000 0.858504 0.568579
        0.429252 0.732633 0.759119
        0.453125 0.694006 0.696619
        0.421875 0.740881 0.618494
        0.500000 0.772131 0.571619
        0.578125 0.740881 0.618494
        0.546875 0.694006 0.696619
        0.616752 0.803381 0.642367
        0.562500 0.880636 0.765198
        0.609375 0.802511 0.728739
        0.569877 0.881506 0.616325
        0.538627 0.928381 0.699658
        0.570748 0.732633 0.759119
        0.500000 0.763013 0.798617
        0.500000 0.856763 0.803825
        0.429252 0.616752 0.687500
        0.616752 0.687500 0.570748
        0.609375 0.609375 0.609375
        0.569877 0.679252 0.500000
        0.500000 0.569877 0.679252
        0.570748 0.616752 0.687500
        0.430123 0.679252 0.500000
        0.500000 0.733504 0.500000
        0.390625 0.609375 0.609375
        0.383248 0.687500 0.570748
        0.858504 0.568579 0.500000
        0.772131 0.571619 0.500000
        0.694006 0.696619 0.453125
        0.740881 0.618494 0.421875
        0.694006 0.696619 0.546875
        0.740881 0.618494 0.578125
        0.880636 0.765198 0.562500
        0.802511 0.728739 0.609375
        0.763013 0.798617 0.500000
        0.856763 0.803825 0.500000
        0.732633 0.759119 0.570748
        0.803381 0.642367 0.616752
        0.928381 0.699658 0.538627
        0.881506 0.616325 0.569877
        0.732633 0.759119 0.429252
        0.880636 0.765198 0.437500
        0.802511 0.728739 0.390625
        0.803381 0.642367 0.383248
        0.881506 0.616325 0.430123
        0.928381 0.699658 0.461373
        0.733504 0.500000 0.500000
        0.687500 0.570748 0.616752
        0.679252 0.500000 0.569877
        0.609375 0.609375 0.390625
        0.616752 0.687500 0.429252
        0.679252 0.500000 0.430123
        0.687500 0.570748 0.383248
        0.568579 0.500000 0.858504
        0.618494 0.421875 0.740881
        0.696619 0.453125 0.694006
        0.571619 0.500000 0.772131
        0.618494 0.578125 0.740881
        0.696619 0.546875 0.694006
        0.765198 0.437500 0.880636
        0.728739 0.390625 0.802511
        0.759119 0.429252 0.732633
        0.699658 0.461373 0.928381
        0.616325 0.430123 0.881506
        0.642367 0.383248 0.803381
        0.765198 0.562500 0.880636
        0.728739 0.609375 0.802511
        0.803825 0.500000 0.856763
        0.798617 0.500000 0.763013
        0.759119 0.570748 0.732633
        0.616325 0.569877 0.881506
        0.699658 0.538627 0.928381
        0.642367 0.616752 0.803381
        0.500000 0.500000 0.733504
        0.500000 0.430123 0.679252
        0.609375 0.390625 0.609375
        0.570748 0.383248 0.687500
        0.687500 0.429252 0.616752
        0.616752 0.929252 0.568579
        0.624129 0.000000 0.616325
        0.577254 0.000000 0.699658
        0.624129 0.000000 0.772575
        0.694006 0.812500 0.618494
        0.694006 0.890625 0.571619
        0.740881 0.812500 0.696619
        0.687500 0.803381 0.759119
        0.687500 0.881506 0.805994
        0.609375 0.928381 0.811202
        0.921004 0.888456 0.570748
        0.881506 0.805994 0.687500
        0.928381 0.811202 0.609375
        0.803381 0.759119 0.687500
        0.732633 0.875871 0.500000
        0.873258 0.881079 0.500000
        0.795133 0.922746 0.500000
        0.929252 0.568579 0.616752
        0.000000 0.616325 0.624129
        0.000000 0.699658 0.577254
        0.000000 0.772575 0.624129
        0.812500 0.618494 0.694006
        0.890625 0.571619 0.694006
        0.812500 0.696619 0.740881
        0.888456 0.570748 0.921004
        0.875871 0.500000 0.732633
        0.922746 0.500000 0.795133
        0.881079 0.500000 0.873258
        0.811202 0.609375 0.928381
        0.805994 0.687500 0.881506
        0.759119 0.687500 0.803381
        0.616325 0.624129 0.000000
        0.699658 0.577254 0.000000
        0.772575 0.624129 0.000000
        0.568579 0.616752 0.929252
        0.618494 0.694006 0.812500
        0.571619 0.694006 0.890625
        0.696619 0.740881 0.812500
        0.500000 0.732633 0.875871
        0.500000 0.795133 0.922746
        0.500000 0.873258 0.881079
        0.570748 0.921004 0.888456
        0.070748 0.767367 0.259119
        0.109375 0.697489 0.228739
        0.062500 0.619364 0.259989
        0.116752 0.696619 0.142367
        0.069877 0.618494 0.126742
        0.038627 0.571619 0.204867
        1.000000 0.641496 0.078996
        0.078125 0.759119 0.118494
        1.000000 0.727869 0.071619
        0.921875 0.759119 0.118494
        0.046875 0.805994 0.196619
        0.953125 0.805994 0.196619
        1.000000 0.736987 0.298617
        1.000000 0.643237 0.298617
        0.929252 0.767367 0.259119
        0.890625 0.697489 0.228739
        0.937500 0.619364 0.259989
        0.883248 0.696619 0.142367
        0.961373 0.571619 0.204867
        0.930123 0.618494 0.126742
        0.070748 0.883248 0.187500
        1.000000 0.930123 0.179252
        0.929252 0.883248 0.187500
        0.109375 0.890625 0.109375
        0.116752 0.812500 0.070748
        0.701383 0.000000 0.356763
        0.701383 1.000000 0.263013
        0.740011 0.937500 0.380636
        0.771261 0.890625 0.302511
        0.740881 0.929252 0.232633
        0.921004 0.000000 0.358504
        0.740011 0.062500 0.380636
        0.771261 0.109375 0.302511
        0.740881 0.070748 0.232633
        0.881506 0.921875 0.240881
        0.881506 0.078125 0.240881
        0.928381 0.000000 0.272131
        0.803381 0.046875 0.194006
        0.803381 0.953125 0.194006
        0.857633 0.116752 0.303381
        0.795133 0.038627 0.428381
        0.873258 0.069877 0.381506
        0.857633 0.883248 0.303381
        0.795133 0.961373 0.428381
        0.873258 0.930123 0.381506
        0.812500 0.070748 0.116752
        0.890625 0.109375 0.109375
        1.000000 0.000000 0.233504
        1.000000 0.069877 0.179252
        0.929252 0.116752 0.187500
        0.500000 0.856763 0.196175
        0.500000 0.763013 0.201383
        0.570748 0.732633 0.240881
        0.437500 0.880636 0.234802
        0.390625 0.802511 0.271261
        0.461373 0.928381 0.300342
        0.430123 0.881506 0.383675
        0.383248 0.803381 0.357633
        0.500000 0.858504 0.431421
        0.429252 0.732633 0.240881
        0.453125 0.694006 0.303381
        0.421875 0.740881 0.381506
        0.578125 0.740881 0.381506
        0.546875 0.694006 0.303381
        0.500000 0.772131 0.428381
        0.562500 0.880636 0.234802
        0.609375 0.802511 0.271261
        0.538627 0.928381 0.300342
        0.569877 0.881506 0.383675
        0.616752 0.803381 0.357633
        0.429252 0.616752 0.312500
        0.500000 0.569877 0.320748
        0.570748 0.616752 0.312500
        0.383248 0.687500 0.429252
        0.390625 0.609375 0.390625
        0.699658 0.538627 0.071619
        0.616325 0.569877 0.118494
        0.642367 0.616752 0.196619
        0.568579 0.500000 0.141496
        0.728739 0.390625 0.197489
        0.765198 0.437500 0.119364
        0.759119 0.429252 0.267367
        0.616325 0.430123 0.118494
        0.699658 0.461373 0.071619
        0.642367 0.383248 0.196619
        0.571619 0.500000 0.227869
        0.618494 0.421875 0.259119
        0.696619 0.453125 0.305994
        0.696619 0.546875 0.305994
        0.618494 0.578125 0.259119
        0.765198 0.562500 0.119364
        0.728739 0.609375 0.197489
        0.759119 0.570748 0.267367
        0.798617 0.500000 0.236987
        0.803825 0.500000 0.143237
        0.609375 0.390625 0.390625
        0.687500 0.429252 0.383248
        0.500000 0.500000 0.266496
        0.500000 0.430123 0.320748
        0.570748 0.383248 0.312500
        0.609375 0.928381 0.188798
        0.687500 0.803381 0.240881
        0.687500 0.881506 0.194006
        0.616752 0.929252 0.431421
        0.577254 0.000000 0.300342
        0.624129 0.000000 0.227425
        0.624129 0.000000 0.383675
        0.694006 0.890625 0.428381
        0.694006 0.812500 0.381506
        0.740881 0.812500 0.303381
        0.921004 0.888456 0.429252
        0.803381 0.759119 0.312500
        0.881506 0.805994 0.312500
        0.928381 0.811202 0.390625
        0.929252 0.568579 0.383248
        0.000000 0.699658 0.422746
        0.000000 0.616325 0.375871
        1.000000 0.772575 0.375871
        0.890625 0.571619 0.305994
        0.812500 0.618494 0.305994
        0.812500 0.696619 0.259119
        0.805994 0.687500 0.118494
        0.811202 0.609375 0.071619
        0.759119 0.687500 0.196619
        0.888456 0.570748 0.078996
        0.922746 0.500000 0.204867
        0.881079 0.500000 0.126742
        0.875871 0.500000 0.267367
        0.618494 0.694006 0.187500
        0.571619 0.694006 0.109375
        0.696619 0.740881 0.187500
        0.568579 0.616752 0.070748
        0.570748 0.921004 0.111544
        0.500000 0.795133 0.077254
        0.500000 0.873258 0.118921
        0.500000 0.732633 0.124129
        0.038627 0.428381 0.795133
        0.069877 0.381506 0.873258
        0.116752 0.303381 0.857633
        0.078125 0.240881 0.881506
        1.000000 0.272131 0.928381
        0.921875 0.240881 0.881506
        0.046875 0.194006 0.803381
        0.953125 0.194006 0.803381
        1.000000 0.356763 0.701383
        0.000000 0.263013 0.701383
        0.937500 0.380636 0.740011
        0.890625 0.302511 0.771261
        0.929252 0.232633 0.740881
        0.883248 0.303381 0.857633
        0.961373 0.428381 0.795133
        0.930123 0.381506 0.873258
        0.000000 0.358504 0.921004
        0.109375 0.302511 0.771261
        0.062500 0.380636 0.740011
        0.070748 0.232633 0.740881
        0.930123 0.179252 1.000000
        0.883248 0.187500 0.929252
        0.000000 0.233504 1.000000
        0.069877 0.179252 1.000000
        0.116752 0.187500 0.929252
        0.109375 0.109375 0.890625
        0.070748 0.116752 0.812500
        0.759119 0.118494 0.921875
        0.805994 0.196619 0.046875
        0.805994 0.196619 0.953125
        0.759119 0.118494 0.078125
        0.727869 0.071619 1.000000
        0.643237 0.298617 1.000000
        0.736987 0.298617 1.000000
        0.767367 0.259119 0.929252
        0.619364 0.259989 0.937500
        0.697489 0.228739 0.890625
        0.571619 0.204867 0.961373
        0.618494 0.126742 0.930123
        0.696619 0.142367 0.883248
        0.641496 0.078996 1.000000
        0.571619 0.204867 0.038627
        0.618494 0.126742 0.069877
        0.696619 0.142367 0.116752
        0.619364 0.259989 0.062500
        0.697489 0.228739 0.109375
        0.767367 0.259119 0.070748
        0.883248 0.187500 0.070748
        0.500000 0.141496 0.568579
        0.421875 0.259119 0.618494
        0.453125 0.305994 0.696619
        0.500000 0.227869 0.571619
        0.578125 0.259119 0.618494
        0.546875 0.305994 0.696619
        0.609375 0.197489 0.728739
        0.562500 0.119364 0.765198
        0.569877 0.118494 0.616325
        0.538627 0.071619 0.699658
        0.616752 0.196619 0.642367
        0.500000 0.236987 0.798617
        0.500000 0.143237 0.803825
        0.570748 0.267367 0.759119
        0.383248 0.196619 0.642367
        0.390625 0.197489 0.728739
        0.437500 0.119364 0.765198
        0.430123 0.118494 0.616325
        0.461373 0.071619 0.699658
        0.429252 0.267367 0.759119
        0.429252 0.383248 0.687500
        0.390625 0.390625 0.609375
        0.569877 0.320748 0.500000
        0.616752 0.312500 0.570748
        0.500000 0.266496 0.500000
        0.430123 0.320748 0.500000
        0.383248 0.312500 0.570748
        0.740881 0.381506 0.421875
        0.694006 0.303381 0.453125
        0.694006 0.303381 0.546875
        0.740881 0.381506 0.578125
        0.772131 0.428381 0.500000
        0.856763 0.196175 0.500000
        0.763013 0.201383 0.500000
        0.880636 0.234802 0.562500
        0.802511 0.271261 0.609375
        0.732633 0.240881 0.570748
        0.858504 0.431421 0.500000
        0.928381 0.300342 0.538627
        0.881506 0.383675 0.569877
        0.803381 0.357633 0.616752
        0.880636 0.234802 0.437500
        0.802511 0.271261 0.390625
        0.732633 0.240881 0.429252
        0.928381 0.300342 0.461373
        0.881506 0.383675 0.430123
        0.803381 0.357633 0.383248
        0.616752 0.312500 0.429252
        0.616752 0.070748 0.568579
        0.694006 0.109375 0.571619
        0.740881 0.187500 0.696619
        0.694006 0.187500 0.618494
        0.687500 0.118494 0.805994
        0.687500 0.196619 0.759119
        0.609375 0.071619 0.811202
        0.921004 0.111544 0.570748
        0.928381 0.188798 0.609375
        0.803381 0.240881 0.687500
        0.881506 0.194006 0.687500
        0.873258 0.118921 0.500000
        0.732633 0.124129 0.500000
        0.795133 0.077254 0.500000
        0.929252 0.431421 0.616752
        0.890625 0.428381 0.694006
        0.812500 0.303381 0.740881
        0.812500 0.381506 0.694006
        0.000000 0.383675 0.624129
        0.000000 0.227425 0.624129
        0.000000 0.300342 0.577254
        0.888456 0.429252 0.921004
        0.811202 0.390625 0.928381
        0.805994 0.312500 0.881506
        0.759119 0.312500 0.803381
        0.568579 0.383248 0.929252
        0.571619 0.305994 0.890625
        0.696619 0.259119 0.812500
        0.618494 0.305994 0.812500
        0.616325 0.375871 0.000000
        0.699658 0.422746 0.000000
        0.772575 0.375871 1.000000
        0.570748 0.078996 0.888456
        0.500000 0.204867 0.922746
        0.500000 0.267367 0.875871
        0.500000 0.126742 0.881079
        0.921875 0.240881 0.118494
        0.078125 0.240881 0.118494
        1.000000 0.272131 0.071619
        0.046875 0.194006 0.196619
        0.953125 0.194006 0.196619
        0.000000 0.358504 0.078996
        1.000000 0.263013 0.298617
        1.000000 0.356763 0.298617
        0.929252 0.232633 0.259119
        0.890625 0.302511 0.228739
        0.937500 0.380636 0.259989
        0.883248 0.303381 0.142367
        0.930123 0.381506 0.126742
        0.961373 0.428381 0.204867
        0.070748 0.232633 0.259119
        0.109375 0.302511 0.228739
        0.062500 0.380636 0.259989
        0.069877 0.381506 0.126742
        0.038627 0.428381 0.204867
        0.116752 0.303381 0.142367
        0.109375 0.109375 0.109375
        0.070748 0.116752 0.187500
        0.116752 0.187500 0.070748
        0.429252 0.267367 0.240881
        0.500000 0.141496 0.431421
        0.421875 0.259119 0.381506
        0.453125 0.305994 0.303381
        0.500000 0.227869 0.428381
        0.546875 0.305994 0.303381
        0.578125 0.259119 0.381506
        0.538627 0.071619 0.300342
        0.616752 0.196619 0.357633
        0.569877 0.118494 0.383675
        0.562500 0.119364 0.234802
        0.609375 0.197489 0.271261
        0.500000 0.143237 0.196175
        0.500000 0.236987 0.201383
        0.570748 0.267367 0.240881
        0.383248 0.196619 0.357633
        0.390625 0.197489 0.271261
        0.437500 0.119364 0.234802
        0.461373 0.071619 0.300342
        0.430123 0.118494 0.383675
        0.429252 0.383248 0.312500
        0.390625 0.390625 0.390625
        0.383248 0.312500 0.429252
        0.616752 0.070748 0.431421
        0.694006 0.187500 0.381506
        0.740881 0.187500 0.303381
        0.694006 0.109375 0.428381
        0.609375 0.071619 0.188798
        0.687500 0.196619 0.240881
        0.687500 0.118494 0.194006
        0.921004 0.111544 0.429252
        0.928381 0.188798 0.390625
        0.881506 0.194006 0.312500
        0.803381 0.240881 0.312500
        0.929252 0.431421 0.383248
        0.812500 0.381506 0.305994
        0.812500 0.303381 0.259119
        0.890625 0.428381 0.305994
        0.000000 0.300342 0.422746
        1.000000 0.227425 0.375871
        0.000000 0.383675 0.375871
        0.888456 0.429252 0.078996
        0.805994 0.312500 0.118494
        0.759119 0.312500 0.196619
        0.811202 0.390625 0.071619
        0.568579 0.383248 0.070748
        0.696619 0.259119 0.187500
        0.571619 0.305994 0.109375
        0.618494 0.305994 0.187500
        0.500000 0.204867 0.077254
        0.500000 0.126742 0.118921
        0.500000 0.267367 0.124129
        0.570748 0.078996 0.111544
        0.263013 0.701383 1.000000
        0.356763 0.701383 0.000000
        0.232633 0.740881 0.929252
        0.302511 0.771261 0.890625
        0.380636 0.740011 0.937500
        0.381506 0.873258 0.930123
        0.428381 0.795133 0.961373
        0.303381 0.857633 0.883248
        0.232633 0.740881 0.070748
        0.302511 0.771261 0.109375
        0.380636 0.740011 0.062500
        0.428381 0.795133 0.038627
        0.381506 0.873258 0.069877
        0.303381 0.857633 0.116752
        0.358504 0.921004 1.000000
        0.194006 0.803381 0.046875
        0.194006 0.803381 0.953125
        0.240881 0.881506 0.078125
        0.272131 0.928381 1.000000
        0.240881 0.881506 0.921875
        0.179252 1.000000 0.069877
        0.187500 0.929252 0.116752
        0.233504 0.000000 1.000000
        0.179252 1.000000 0.930123
        0.187500 0.929252 0.883248
        0.259119 0.929252 0.767367
        0.228739 0.890625 0.697489
        0.259989 0.937500 0.619364
        0.298617 1.000000 0.736987
        0.298617 1.000000 0.643237
        0.142367 0.883248 0.696619
        0.204867 0.961373 0.571619
        0.126742 0.930123 0.618494
        0.078996 1.000000 0.641496
        0.259119 0.070748 0.767367
        0.228739 0.109375 0.697489
        0.259989 0.062500 0.619364
        0.142367 0.116752 0.696619
        0.126742 0.069877 0.618494
        0.204867 0.038627 0.571619
        0.196619 0.046875 0.805994
        0.196619 0.953125 0.805994
        0.118494 0.078125 0.759119
        0.071619 1.000000 0.727869
        0.118494 0.921875 0.759119
        0.187500 0.070748 0.883248
        0.141496 0.568579 0.500000
        0.267367 0.759119 0.570748
        0.197489 0.728739 0.609375
        0.119364 0.765198 0.562500
        0.236987 0.798617 0.500000
        0.143237 0.803825 0.500000
        0.118494 0.616325 0.569877
        0.071619 0.699658 0.538627
        0.196619 0.642367 0.616752
        0.197489 0.728739 0.390625
        0.119364 0.765198 0.437500
        0.267367 0.759119 0.429252
        0.196619 0.642367 0.383248
        0.118494 0.616325 0.430123
        0.071619 0.699658 0.461373
        0.227869 0.571619 0.500000
        0.259119 0.618494 0.421875
        0.305994 0.696619 0.453125
        0.259119 0.618494 0.578125
        0.305994 0.696619 0.546875
        0.312500 0.570748 0.616752
        0.320748 0.500000 0.569877
        0.312500 0.570748 0.383248
        0.320748 0.500000 0.430123
        0.266496 0.500000 0.500000
        0.303381 0.453125 0.694006
        0.381506 0.421875 0.740881
        0.381506 0.578125 0.740881
        0.303381 0.546875 0.694006
        0.428381 0.500000 0.772131
        0.196175 0.500000 0.856763
        0.201383 0.500000 0.763013
        0.240881 0.570748 0.732633
        0.234802 0.562500 0.880636
        0.271261 0.609375 0.802511
        0.300342 0.538627 0.928381
        0.357633 0.616752 0.803381
        0.383675 0.569877 0.881506
        0.431421 0.500000 0.858504
        0.240881 0.429252 0.732633
        0.234802 0.437500 0.880636
        0.271261 0.390625 0.802511
        0.300342 0.461373 0.928381
        0.357633 0.383248 0.803381
        0.383675 0.430123 0.881506
        0.312500 0.429252 0.616752
        0.259119 0.812500 0.696619
        0.305994 0.812500 0.618494
        0.305994 0.890625 0.571619
        0.383248 0.929252 0.568579
        0.312500 0.881506 0.805994
        0.312500 0.803381 0.759119
        0.390625 0.928381 0.811202
        0.422746 0.000000 0.699658
        0.375871 0.000000 0.616325
        0.375871 0.000000 0.772575
        0.196619 0.759119 0.687500
        0.071619 0.811202 0.609375
        0.118494 0.805994 0.687500
        0.078996 0.888456 0.570748
        0.204867 0.922746 0.500000
        0.126742 0.881079 0.500000
        0.267367 0.875871 0.500000
        0.070748 0.568579 0.616752
        0.187500 0.696619 0.740881
        0.187500 0.618494 0.694006
        0.109375 0.571619 0.694006
        0.194006 0.687500 0.881506
        0.240881 0.687500 0.803381
        0.188798 0.609375 0.928381
        0.111544 0.570748 0.921004
        0.077254 0.500000 0.795133
        0.124129 0.500000 0.732633
        0.118921 0.500000 0.873258
        0.303381 0.740881 0.812500
        0.381506 0.694006 0.812500
        0.428381 0.694006 0.890625
        0.431421 0.616752 0.929252
        0.227425 0.624129 0.000000
        0.300342 0.577254 0.000000
        0.383675 0.624129 0.000000
        0.429252 0.921004 0.888456
        0.118494 0.921875 0.240881
        0.196619 0.046875 0.194006
        0.196619 0.953125 0.194006
        0.118494 0.078125 0.240881
        0.071619 0.000000 0.272131
        0.298617 1.000000 0.263013
        0.298617 1.000000 0.356763
        0.259119 0.929252 0.232633
        0.228739 0.890625 0.302511
        0.259989 0.937500 0.380636
        0.142367 0.883248 0.303381
        0.126742 0.930123 0.381506
        0.204867 0.961373 0.428381
        0.078996 0.000000 0.358504
        0.259119 0.070748 0.232633
        0.228739 0.109375 0.302511
        0.259989 0.062500 0.380636
        0.142367 0.116752 0.303381
        0.126742 0.069877 0.381506
        0.204867 0.038627 0.428381
        0.187500 0.070748 0.116752
        0.357633 0.383248 0.196619
        0.383675 0.430123 0.118494
        0.300342 0.461373 0.071619
        0.431421 0.500000 0.141496
        0.303381 0.453125 0.305994
        0.381506 0.421875 0.259119
        0.303381 0.546875 0.305994
        0.381506 0.578125 0.259119
        0.428381 0.500000 0.227869
        0.240881 0.570748 0.267367
        0.234802 0.562500 0.119364
        0.271261 0.609375 0.197489
        0.196175 0.500000 0.143237
        0.201383 0.500000 0.236987
        0.300342 0.538627 0.071619
        0.357633 0.616752 0.196619
        0.383675 0.569877 0.118494
        0.234802 0.437500 0.119364
        0.271261 0.390625 0.197489
        0.240881 0.429252 0.267367
        0.312500 0.429252 0.383248
        0.383248 0.929252 0.431421
        0.305994 0.890625 0.428381
        0.259119 0.812500 0.303381
        0.305994 0.812500 0.381506
        0.312500 0.803381 0.240881
        0.312500 0.881506 0.194006
        0.390625 0.928381 0.188798
        0.375871 0.000000 0.227425
        0.422746 0.000000 0.300342
        0.375871 0.000000 0.383675
        0.078996 0.888456 0.429252
        0.196619 0.759119 0.312500
        0.118494 0.805994 0.312500
        0.071619 0.811202 0.390625
        0.070748 0.568579 0.383248
        0.109375 0.571619 0.305994
        0.187500 0.696619 0.259119
        0.187500 0.618494 0.305994
        0.111544 0.570748 0.078996
        0.188798 0.609375 0.071619
        0.240881 0.687500 0.196619
        0.194006 0.687500 0.118494
        0.118921 0.500000 0.126742
        0.124129 0.500000 0.267367
        0.077254 0.500000 0.204867
        0.431421 0.616752 0.070748
        0.303381 0.740881 0.187500
        0.428381 0.694006 0.109375
        0.381506 0.694006 0.187500
        0.429252 0.921004 0.111544
        0.303381 0.142367 0.116752
        0.428381 0.204867 0.038627
        0.381506 0.126742 0.069877
        0.302511 0.228739 0.109375
        0.380636 0.259989 0.062500
        0.232633 0.259119 0.070748
        0.356763 0.298617 1.000000
        0.263013 0.298617 1.000000
        0.380636 0.259989 0.937500
        0.302511 0.228739 0.890625
        0.232633 0.259119 0.929252
        0.358504 0.078996 1.000000
        0.240881 0.118494 0.921875
        0.240881 0.118494 0.078125
        0.272131 0.071619 1.000000
        0.194006 0.196619 0.046875
        0.194006 0.196619 0.953125
        0.428381 0.204867 0.961373
        0.381506 0.126742 0.930123
        0.303381 0.142367 0.883248
        0.141496 0.431421 0.500000
        0.259119 0.381506 0.421875
        0.305994 0.303381 0.453125
        0.305994 0.303381 0.546875
        0.259119 0.381506 0.578125
        0.227869 0.428381 0.500000
        0.197489 0.271261 0.609375
        0.119364 0.234802 0.562500
        0.267367 0.240881 0.570748
        0.236987 0.201383 0.500000
        0.143237 0.196175 0.500000
        0.196619 0.357633 0.616752
        0.118494 0.383675 0.569877
        0.071619 0.300342 0.538627
        0.267367 0.240881 0.429252
        0.197489 0.271261 0.390625
        0.119364 0.234802 0.437500
        0.196619 0.357633 0.383248
        0.071619 0.300342 0.461373
        0.118494 0.383675 0.430123
        0.383248 0.070748 0.568579
        0.259119 0.187500 0.696619
        0.305994 0.187500 0.618494
        0.305994 0.109375 0.571619
        0.312500 0.196619 0.759119
        0.312500 0.118494 0.805994
        0.390625 0.071619 0.811202
        0.078996 0.111544 0.570748
        0.196619 0.240881 0.687500
        0.071619 0.188798 0.609375
        0.118494 0.194006 0.687500
        0.126742 0.118921 0.500000
        0.267367 0.124129 0.500000
        0.204867 0.077254 0.500000
        0.187500 0.303381 0.740881
        0.187500 0.381506 0.694006
        0.109375 0.428381 0.694006
        0.070748 0.431421 0.616752
        0.111544 0.429252 0.921004
        0.240881 0.312500 0.803381
        0.188798 0.390625 0.928381
        0.194006 0.312500 0.881506
        0.300342 0.422746 0.000000
        0.227425 0.375871 1.000000
        0.383675 0.375871 0.000000
        0.431421 0.383248 0.929252
        0.428381 0.305994 0.890625
        0.303381 0.259119 0.812500
        0.381506 0.305994 0.812500
        0.429252 0.078996 0.888456
        0.312500 0.196619 0.240881
        0.312500 0.118494 0.194006
        0.390625 0.071619 0.188798
        0.383248 0.070748 0.431421
        0.259119 0.187500 0.303381
        0.305994 0.109375 0.428381
        0.305994 0.187500 0.381506
        0.078996 0.111544 0.429252
        0.196619 0.240881 0.312500
        0.071619 0.188798 0.390625
        0.118494 0.194006 0.312500
        0.070748 0.431421 0.383248
        0.187500 0.303381 0.259119
        0.109375 0.428381 0.305994
        0.187500 0.381506 0.305994
        0.240881 0.312500 0.196619
        0.194006 0.312500 0.118494
        0.188798 0.390625 0.071619
        0.111544 0.429252 0.078996
        0.303381 0.259119 0.187500
        0.381506 0.305994 0.187500
        0.428381 0.305994 0.109375
        0.431421 0.383248 0.070748
        0.429252 0.078996 0.111544
        0.500000 0.421875 0.046875
        0.500000 0.500000 0.093750
        0.500000 0.578125 0.046875
        0.500000 0.421875 0.953125
        0.500000 0.578125 0.953125
        0.500000 0.500000 0.906250
        0.093750 0.500000 0.500000
        0.046875 0.500000 0.421875
        0.046875 0.500000 0.578125
        0.953125 0.500000 0.421875
        0.953125 0.500000 0.578125
        0.906250 0.500000 0.500000
        0.500000 0.093750 0.500000
        0.421875 0.046875 0.500000
        0.578125 0.046875 0.500000
        0.421875 0.953125 0.500000
        0.578125 0.953125 0.500000
        0.500000 0.906250 0.500000
        0.927083 0.453125 1.000000
        0.927083 0.546875 1.000000
        1.000000 0.927083 0.453125
        0.000000 0.927083 0.546875
        0.453125 1.000000 0.927083
        0.546875 1.000000 0.927083
        0.072917 0.546875 1.000000
        0.072917 0.453125 1.000000
        0.000000 0.072917 0.546875
        1.000000 0.072917 0.453125
        0.453125 1.000000 0.072917
        0.546875 1.000000 0.072917
        1.000000 0.000000 0.406250
        0.000000 1.000000 0.593750
        0.406250 1.000000 0.000000
        0.593750 1.000000 0.000000
        0.000000 0.406250 1.000000
        0.000000 0.593750 1.000000
        """

        self.cell = cellvectors(a=2.75276384094,
                           b=2.75276384094,
                           c=2.75276384094)
