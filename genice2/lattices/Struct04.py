# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (16,0,0,8,)
"""
from genice2.cell import cellvectors
import genice2.lattices
desc = {"ref": {"SpaceFullerene": 'Sikiric 2010'},
        "usage": "No options available.",
        "brief": "A space fullerene."
        }


class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.pairs = """
        37 101
        72 133
        121 134
        99 56
        84 6
        124 101
        104 133
        9 118
        119 79
        29 67
        2 106
        83 113
        18 50
        82 116
        7 78
        9 25
        117 84
        11 9
        98 79
        68 2
        130 77
        67 3
        132 79
        36 63
        105 91
        45 101
        135 92
        19 119
        40 59
        90 47
        27 8
        36 81
        83 99
        37 58
        74 99
        12 80
        60 48
        25 32
        18 131
        103 95
        65 91
        82 15
        106 44
        110 48
        97 26
        53 124
        80 64
        127 102
        45 54
        59 129
        7 41
        34 15
        111 93
        70 135
        107 100
        16 54
        61 108
        0 64
        94 86
        7 58
        109 125
        68 117
        66 15
        112 52
        90 85
        120 48
        53 114
        23 56
        2 6
        22 71
        123 33
        8 5
        92 128
        68 24
        65 95
        117 5
        135 25
        50 66
        45 72
        46 73
        10 42
        44 10
        86 115
        2 36
        65 133
        120 28
        82 71
        110 25
        1 69
        11 111
        14 39
        77 78
        57 134
        62 116
        36 30
        63 35
        112 115
        51 95
        14 90
        88 64
        122 10
        12 131
        35 20
        90 115
        107 76
        72 51
        112 65
        30 83
        59 55
        94 103
        67 123
        23 30
        119 103
        49 33
        8 63
        16 78
        7 0
        21 3
        59 58
        24 23
        121 124
        18 22
        83 97
        128 32
        132 109
        22 40
        0 129
        31 80
        80 96
        44 62
        121 75
        85 28
        122 26
        89 69
        16 134
        111 98
        46 60
        44 63
        24 126
        94 54
        29 43
        40 114
        93 107
        130 95
        119 73
        74 128
        29 89
        110 84
        123 62
        132 52
        24 122
        96 55
        39 32
        89 129
        27 21
        79 28
        19 118
        38 96
        106 33
        51 47
        109 118
        64 71
        31 78
        76 98
        11 17
        133 86
        70 81
        38 57
        52 103
        56 5
        84 99
        61 129
        1 22
        77 94
        113 92
        46 98
        120 76
        34 20
        60 125
        93 32
        121 130
        47 73
        17 127
        57 41
        17 135
        87 3
        53 58
        39 118
        13 35
        102 56
        81 5
        38 105
        53 91
        87 42
        122 20
        61 88
        37 16
        41 91
        49 8
        13 26
        61 66
        72 109
        13 21
        45 52
        21 33
        13 30
        97 6
        100 113
        34 108
        131 88
        50 69
        76 92
        100 126
        60 9
        123 15
        93 127
        49 68
        50 67
        75 115
        17 120
        31 114
        4 108
        97 126
        75 54
        106 26
        70 102
        77 104
        40 12
        111 14
        125 47
        11 85
        104 101
        31 105
        29 108
        74 48
        12 89
        18 87
        35 4
        4 116
        70 6
        124 112
        14 132
        74 127
        125 28
        117 100
        130 41
        75 51
        46 128
        131 43
        113 81
        66 55
        27 10
        55 71
        23 27
        69 96
        34 42
        102 126
        38 37
        1 0
        43 42
        4 3
        39 73
        19 85
        19 86
        134 114
        82 87
        110 107
        88 116
        104 105
        49 20
        57 1
        62 43
        """

        self.waters = """
        0.03915 0.71085 0.84348
        0.28915 0.96085 0.84348
        0.75 0.875 0.5
        0.28915 0.96085 0.65653
        0.03915 0.71085 0.65653
        0.66667 0.33334 0.46828
        0.58334 0.79167 0.4377
        0.0 0.625 0.90598
        0.66667 0.33334 0.53907
        0.58837 0.79418 0.25
        0.375 0.375 0.59403
        0.45582 0.54418 0.25
        0.46085 0.53916 0.80421
        0.20833 0.79167 0.55496
        0.12249 0.24497 0.1873
        0.87752 0.12249 0.6873
        0.66667 0.33333 0.96094
        0.53916 0.46085 0.30421
        0.41164 0.20582 0.75
        0.71085 0.6717 0.15653
        0.0 0.375 0.59403
        0.375 0.0 0.59403
        0.24497 0.12249 0.8127
        0.25 0.125 0.5
        0.125 0.25 0.5
        0.75503 0.87752 0.3127
        0.33333 0.66667 0.53173
        0.41667 0.20833 0.5623
        0.66667 0.33333 0.21923
        0.33334 0.66667 0.71923
        0.125 0.875 0.5
        0.625 0.625 0.90598
        0.0 0.0 0.28927
        0.625 0.0 0.59403
        0.03916 0.32831 0.65653
        0.0 0.625 0.59403
        0.875 0.75 0.5
        0.79167 0.20833 0.9377
        0.625 0.0 0.90598
        0.0 0.0 0.21074
        0.28915 0.32831 0.84348
        0.20834 0.79167 0.94504
        0.28915 0.32831 0.65653
        0.46085 0.53915 0.69579
        0.625 0.625 0.59403
        0.79167 0.20833 0.05496
        0.20582 0.79418 0.25
        0.32831 0.03915 0.15653
        0.53915 0.07831 0.30421
        0.79167 0.20833 0.5623
        0.54418 0.08836 0.75
        0.375 0.0 0.09403
        0.0 0.375 0.09403
        0.20833 0.41667 0.94504
        0.66667 0.33333 0.03173
        0.87752 0.12249 0.8127
        0.41667 0.20834 0.44504
        0.375 0.0 0.90598
        0.0 0.375 0.90598
        0.03915 0.32831 0.84348
        0.45582 0.91164 0.25
        0.91164 0.45582 0.75
        0.67169 0.71085 0.65653
        0.79167 0.58334 0.5623
        0.87752 0.75503 0.8127
        0.33334 0.66667 0.03907
        0.79418 0.20582 0.75
        0.46085 0.9217 0.69579
        0.875 0.125 0.5
        0.46085 0.9217 0.80421
        0.625 0.625 0.40598
        0.0 0.0 0.78927
        0.625 0.0 0.09403
        0.12249 0.87752 0.1873
        0.32831 0.03915 0.34348
        0.41667 0.20833 0.05496
        0.9217 0.46085 0.30421
        0.875 0.75 0.0
        0.79167 0.58334 0.9377
        0.92169 0.46085 0.19579
        0.6717 0.71085 0.84348
        0.79167 0.58334 0.44504
        0.0 0.0 0.71074
        0.20833 0.79167 0.4377
        0.625 0.0 0.40598
        0.53915 0.46085 0.19579
        0.625 0.625 0.09403
        0.24497 0.12249 0.6873
        0.79418 0.58837 0.75
        0.33333 0.66667 0.78078
        0.3283 0.28915 0.15653
        0.33334 0.66667 0.96828
        0.96085 0.6717 0.34348
        0.12249 0.24497 0.3127
        0.79167 0.58333 0.05496
        0.20834 0.79167 0.0623
        0.6717 0.96085 0.84348
        0.33334 0.66667 0.46094
        0.08836 0.54418 0.25
        0.375 0.0 0.40598
        0.0 0.375 0.40598
        0.875 0.125 0.0
        0.375 0.375 0.40598
        0.0 0.625 0.09403
        0.75 0.875 0.0
        0.58334 0.79167 0.94504
        0.58333 0.79167 0.55496
        0.96085 0.28915 0.34348
        0.07831 0.53915 0.69579
        0.71085 0.03915 0.15653
        0.71085 0.03915 0.34348
        0.20582 0.41164 0.25
        0.20833 0.41667 0.0623
        0.0 0.625 0.40598
        0.375 0.375 0.90598
        0.375 0.375 0.09403
        0.87752 0.75503 0.6873
        0.79167 0.20833 0.44504
        0.75503 0.87752 0.1873
        0.96085 0.6717 0.15653
        0.66667 0.33333 0.28078
        0.25 0.125 0.0
        0.20833 0.41667 0.55496
        0.6717 0.96085 0.65653
        0.125 0.25 0.0
        0.53915 0.07831 0.19579
        0.20834 0.41667 0.4377
        0.32831 0.28915 0.34348
        0.12249 0.87752 0.3127
        0.07831 0.53915 0.80421
        0.125 0.875 0.0
        0.54418 0.45582 0.75
        0.96085 0.28916 0.15653
        0.58334 0.79167 0.0623
        0.41667 0.20833 0.9377
        0.71085 0.6717 0.34348
        """

        self.coord = "relative"

        self.cages = """
        12 0.5 0.0 0.0
        12 -0.15661 -0.31322 0.25
        12 0.33333 0.66667 -0.1269
        12 0.15661 0.31322 -0.25
        12 -0.66667 -0.33333 0.6269
        16 0.33333 0.66667 0.15626
        12 0.31322 0.15661 0.25
        16 0.0 0.0 0.59294
        12 0.5 0.5 0.5
        12 -0.5 0.0 -0.5
        16 0.0 0.0 -0.09294
        12 -0.15661 0.15661 0.25
        12 0.0 0.5 0.0
        12 0.0 -0.5 0.5
        16 0.0 0.0 0.09294
        16 0.0 0.0 -0.59294
        12 0.15661 -0.15661 0.75
        16 0.66667 0.33333 -0.15626
        12 0.66667 0.33333 0.1269
        12 -0.31322 -0.15661 0.75
        12 -0.5 -0.5 0.0
        16 -0.33333 -0.66667 0.65626
        16 -0.66667 -0.33333 -0.65626
        12 -0.33333 -0.66667 -0.6269
        """

        self.bondlen = 3

        self.cell = """
        13.698097440028898 0.0 0.0
        -6.849048720014445 11.862900366579614 0.0
        2.7457556179095828e-15 4.755788235387074e-15 44.84159220146225
        """

        self.density = 0.5578770598708237

        self.cell = cellvectors(a=13.698097440028898,
                                b=13.6980974400289,
                                c=44.84159220146225,
                                C=119.99999999999999)
