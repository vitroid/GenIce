# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (10,4,4,2,)
"""
desc={"ref": {"SpaceFullerene": 'Sikiric 2010'},
      "usage": "No options available.",
      "brief": "A space fullerene."
      }

import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.pairs="""
        70 45
        110 8
        4 68
        62 50
        24 89
        95 87
        105 83
        97 76
        11 64
        99 78
        0 93
        34 69
        111 42
        37 62
        48 25
        71 44
        46 36
        97 38
        37 101
        23 20
        4 47
        30 22
        63 56
        73 81
        41 36
        4 5
        84 87
        86 108
        109 31
        24 96
        0 14
        9 50
        53 107
        1 33
        92 109
        28 92
        35 82
        4 39
        12 107
        62 54
        46 6
        30 100
        92 97
        99 94
        65 90
        66 89
        74 39
        17 26
        57 95
        63 20
        9 77
        65 113
        63 27
        58 77
        43 66
        111 63
        75 60
        81 32
        14 51
        97 22
        103 100
        8 47
        66 84
        73 54
        32 69
        60 31
        84 51
        108 20
        28 43
        104 87
        82 6
        55 38
        16 58
        74 2
        14 69
        104 42
        79 21
        12 59
        96 44
        101 34
        25 78
        19 107
        12 57
        25 54
        31 41
        32 42
        95 78
        79 99
        15 16
        13 46
        72 108
        84 50
        12 29
        31 6
        18 57
        98 33
        15 9
        30 73
        70 71
        90 91
        8 7
        60 101
        22 34
        106 68
        93 47
        43 9
        2 100
        11 93
        62 105
        44 112
        103 67
        18 33
        42 54
        52 112
        27 30
        106 113
        35 88
        70 52
        40 102
        48 70
        101 73
        17 110
        58 55
        112 20
        0 90
        1 91
        74 106
        110 41
        61 83
        78 50
        14 104
        15 105
        37 61
        41 68
        24 80
        28 6
        60 64
        75 68
        56 22
        1 39
        19 79
        89 16
        65 112
        82 64
        79 77
        33 38
        3 8
        113 85
        36 26
        100 76
        28 88
        23 96
        16 87
        111 59
        111 95
        65 51
        10 110
        53 99
        13 71
        21 76
        98 94
        82 102
        80 45
        94 57
        74 21
        109 34
        48 53
        85 2
        72 67
        56 32
        17 94
        49 55
        75 93
        72 2
        75 61
        53 26
        98 21
        23 45
        90 108
        61 113
        37 51
        5 36
        52 85
        25 27
        29 7
        44 66
        10 106
        86 83
        103 81
        86 85
        35 103
        109 3
        11 69
        40 46
        107 5
        48 59
        49 102
        52 27
        67 38
        13 43
        92 49
        40 7
        17 29
        49 89
        7 5
        24 40
        56 67
        96 29
        11 3
        104 105
        86 45
        1 72
        18 58
        35 55
        76 88
        77 88
        0 83
        91 47
        10 91
        19 18
        19 39
        3 102
        71 26
        13 80
        15 80
        10 98
        23 59
        64 81
        """

        self.waters="""
        0.00328 0.69394 0.8029
        0.66618 0.69397 0.58101
        0.63813 0.375 0.69928
        0.07507 0.81897 0.52457
        0.80606 0.5 0.50388
        0.80206 0.5 0.37005
        0.11959 0.30603 0.42749
        0.87785 0.69397 0.35055
        0.92494 0.81897 0.47544
        0.23687 0.30606 0.17572
        0.7592 0.0 0.5475
        0.12216 0.69397 0.64945
        0.64315 0.69397 0.19548
        0.99672 0.30606 0.1971
        0.13418 0.81894 0.85532
        0.14692 0.5 0.10328
        0.23687 0.69394 0.17572
        0.7635 0.0 0.30138
        0.493 0.69397 0.32917
        0.57798 0.5 0.35716
        0.71985 0.81894 0.91584
        0.54994 0.18103 0.47544
        0.39727 0.0 0.70712
        0.762 0.69394 0.03246
        0.99672 0.69394 0.1971
        0.5 0.25 0.0
        0.79076 0.18103 0.24512
        0.55694 0.125 0.89627
        0.20954 0.18103 0.37493
        0.79076 0.81897 0.24512
        0.46909 0.18103 0.77627
        0.07507 0.18103 0.52457
        0.35685 0.69397 0.80453
        0.54994 0.81897 0.47544
        0.2365 0.0 0.69863
        0.31894 0.5 0.49612
        0.87785 0.30603 0.35055
        0.13418 0.18106 0.85532
        0.45007 0.81897 0.52457
        0.68106 0.5 0.50388
        0.99672 0.625 0.3221
        0.92494 0.18103 0.47544
        0.3849 0.625 0.93625
        0.13448 0.18106 0.22536
        0.91363 0.0 0.10051
        0.84908 0.5 0.01289
        0.99672 0.375 0.3221
        0.88041 0.69397 0.57251
        0.61511 0.375 0.06375
        0.20954 0.81897 0.37493
        0.28015 0.18106 0.08416
        0.08637 0.0 0.8995
        0.71985 0.18106 0.91584
        0.64315 0.30603 0.19548
        0.3849 0.375 0.93625
        0.33383 0.69397 0.41899
        0.46909 0.81897 0.77627
        0.53091 0.81897 0.22374
        0.36187 0.625 0.30072
        0.61511 0.625 0.06375
        0.12216 0.30603 0.64945
        0.00328 0.30606 0.8029
        0.23801 0.30606 0.96755
        0.55694 0.875 0.89627
        0.19794 0.5 0.62995
        0.90933 0.0 0.84663
        0.09067 0.0 0.15337
        0.507 0.69397 0.67084
        0.88041 0.30603 0.57251
        0.20924 0.81897 0.75488
        0.762 0.30606 0.03246
        0.86582 0.18106 0.14468
        0.63813 0.625 0.69928
        0.35685 0.30603 0.80453
        0.66618 0.30603 0.58101
        0.00328 0.375 0.6779
        0.45007 0.18103 0.52457
        0.36187 0.375 0.30072
        0.44306 0.125 0.10373
        0.493 0.30603 0.32917
        0.996 0.5 0.11617
        0.34483 0.5 0.72366
        0.19394 0.5 0.49612
        0.004 0.5 0.88383
        0.17704 0.0 0.05287
        0.76313 0.30606 0.82428
        0.85308 0.5 0.89672
        0.28015 0.81894 0.08416
        0.33383 0.30603 0.41899
        0.13448 0.81894 0.22536
        0.86552 0.81894 0.77464
        0.79046 0.81897 0.62507
        0.2408 0.0 0.4525
        0.00328 0.625 0.6779
        0.60274 0.0 0.29288
        0.44306 0.875 0.10373
        0.86582 0.81894 0.14468
        0.40013 0.0 0.54913
        0.59987 0.0 0.45087
        0.53091 0.18103 0.22374
        0.507 0.30603 0.67084
        0.20924 0.18103 0.75488
        0.11959 0.69397 0.42749
        0.42202 0.5 0.64284
        0.23801 0.69394 0.96755
        0.15092 0.5 0.98712
        0.79046 0.18103 0.62507
        0.65517 0.5 0.27634
        0.76313 0.69394 0.82428
        0.15013 0.0 0.54913
        0.84987 0.0 0.45087
        0.5 0.75 0.0
        0.82296 0.0 0.94713
        0.86552 0.18106 0.77464
        """

        self.coord= "relative"

        self.cages="""
        12 0.80877 0.5 0.17671
        14 0.0 -0.22423 0.0
        16 -0.35981 0.0 -0.28975
        12 0.30026 0.22411 0.59826
        15 -0.00287 0.0 -0.32374
        12 -0.80877 0.5 -0.17671
        12 -0.17523 0.5 -0.28796
        15 -0.41244 0.5 -0.12517
        12 -0.30026 -0.22411 -0.59826
        12 0.5 0.5 0.5
        15 0.00287 0.0 0.32374
        14 0.65165 0.0 0.07828
        12 0.30026 -0.22411 0.59826
        12 0.0 0.5 0.5
        14 0.0 0.22423 0.0
        14 -0.65165 0.0 -0.07828
        12 0.17523 0.5 0.28796
        16 0.35981 0.0 0.28975
        12 -0.30026 0.22411 -0.59826
        15 0.41244 0.5 0.12517
        """

        self.bondlen = 3


        self.cell = """
        22.080000000000002 0.0 0.0
        8.99051812062174e-16 14.682630333711383 0.0
        -8.26777972121982 1.9962906425462054e-15 22.886529726450306
        """

        self.density = 0.4592548228777292



        self.cell = cellvectors(a=22.080000000000002,
                           b=14.682630333711383,
                           c=24.334120580746433,
                           B=109.86230000000002)
