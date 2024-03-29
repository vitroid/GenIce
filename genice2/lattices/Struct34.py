# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (12,8,8,0,)
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
        129 46
        15 96
        109 30
        4 151
        19 128
        87 99
        32 100
        20 154
        65 36
        72 6
        113 35
        132 156
        20 44
        135 37
        121 72
        122 32
        35 44
        38 28
        111 10
        115 43
        95 108
        86 21
        99 94
        152 10
        50 96
        156 137
        9 63
        122 159
        101 57
        15 130
        124 128
        117 131
        76 158
        107 102
        104 18
        23 16
        129 134
        121 110
        104 6
        93 21
        2 116
        12 138
        76 138
        89 154
        69 82
        95 22
        105 148
        9 53
        105 149
        20 34
        1 82
        143 100
        3 6
        158 33
        127 133
        78 128
        142 43
        14 70
        108 71
        0 36
        124 79
        124 80
        77 119
        136 59
        89 117
        88 118
        22 73
        112 155
        142 53
        118 151
        137 120
        103 48
        92 68
        116 38
        105 133
        116 53
        51 66
        26 81
        121 49
        40 141
        136 55
        142 120
        33 84
        43 80
        157 25
        93 120
        3 125
        145 158
        74 125
        106 42
        99 12
        97 149
        98 155
        69 36
        48 24
        91 72
        156 21
        8 84
        5 49
        13 114
        74 62
        92 134
        64 70
        13 95
        157 149
        87 9
        148 128
        67 102
        25 0
        109 58
        132 27
        115 47
        87 2
        134 100
        110 15
        103 100
        95 34
        41 84
        11 54
        4 123
        78 57
        37 25
        17 57
        158 56
        50 66
        67 117
        92 11
        35 123
        110 26
        144 26
        17 143
        144 18
        152 31
        98 96
        41 148
        12 141
        24 90
        142 153
        131 85
        87 40
        88 42
        69 59
        61 58
        50 81
        126 4
        19 101
        112 139
        46 101
        85 139
        121 86
        107 30
        97 41
        8 115
        113 77
        76 3
        77 1
        113 61
        109 89
        111 37
        131 107
        19 115
        61 154
        23 137
        147 151
        78 141
        73 59
        60 73
        31 90
        30 91
        92 34
        144 66
        114 82
        97 60
        54 64
        145 38
        127 141
        37 65
        13 31
        7 108
        29 159
        61 136
        137 62
        45 41
        33 43
        76 80
        44 117
        126 29
        130 66
        39 46
        154 11
        109 139
        104 80
        129 14
        31 71
        102 51
        98 4
        112 106
        2 148
        20 55
        105 40
        112 81
        150 52
        120 104
        29 140
        27 63
        13 48
        44 119
        70 59
        123 58
        54 122
        130 86
        113 69
        67 159
        79 18
        38 16
        74 16
        114 65
        93 94
        145 63
        143 24
        99 146
        39 24
        156 107
        88 136
        26 28
        78 45
        152 52
        103 119
        135 60
        94 18
        111 149
        111 150
        3 21
        71 101
        147 42
        77 11
        39 133
        145 125
        74 139
        42 140
        106 58
        5 91
        93 138
        84 52
        22 55
        23 53
        32 0
        79 28
        54 140
        1 48
        16 81
        7 83
        7 82
        108 68
        8 83
        131 123
        132 125
        144 72
        19 133
        151 159
        119 122
        118 32
        2 79
        68 46
        64 134
        102 86
        7 73
        146 28
        1 68
        97 83
        98 85
        23 27
        147 49
        75 36
        126 106
        75 118
        75 55
        75 35
        103 34
        127 90
        62 6
        10 83
        143 25
        146 63
        14 60
        126 50
        135 45
        150 40
        88 64
        5 140
        62 30
        132 85
        9 56
        15 146
        56 52
        67 49
        110 155
        57 90
        91 51
        17 129
        153 138
        70 0
        157 14
        22 65
        130 94
        157 39
        45 150
        127 47
        89 5
        29 51
        47 152
        33 116
        47 153
        8 71
        114 10
        135 17
        56 153
        147 155
        27 96
        12 124
        """

        self.waters = """
        0.375 0.09721 0.35098
        0.5 0.05554 0.875
        0.625 0.70096 0.36774
        0.0 0.55554 0.875
        0.3125 0.32598 0.46879
        0.8125 0.32598 0.03122
        0.8125 0.52777 0.9375
        0.6875 0.98528 0.6875
        0.75 0.84808 0.75
        0.3125 0.67402 0.53122
        0.375 0.90279 0.64902
        0.6875 0.17402 0.96879
        0.125 0.70096 0.13227
        0.1875 0.98528 0.8125
        0.6875 0.01473 0.3125
        0.1875 0.51473 0.3125
        0.6875 0.52777 0.5625
        0.0 0.96571 0.21872
        0.625 0.59721 0.14902
        0.6875 0.82598 0.03122
        0.0 0.18054 0.79804
        0.1875 0.52777 0.9375
        0.0 0.05554 0.625
        0.5 0.55554 0.625
        0.3125 0.97223 0.0625
        0.3125 0.01473 0.3125
        0.8125 0.51473 0.3125
        0.3125 0.52777 0.5625
        0.875 0.59721 0.35098
        0.5 0.31946 0.20196
        0.625 0.40279 0.85098
        0.125 0.90279 0.85098
        0.25 0.15193 0.25
        0.8125 0.72306 0.66348
        0.0 0.11026 0.89902
        0.3125 0.22306 0.66348
        0.3125 0.10373 0.5072
        0.1875 0.97223 0.4375
        0.8125 0.60373 0.5072
        0.5 0.94446 0.125
        0.3125 0.77695 0.33652
        0.8125 0.82598 0.46879
        0.8125 0.27695 0.33652
        0.6875 0.72306 0.83652
        0.1875 0.22306 0.83652
        0.0 0.85293 0.38227
        0.6875 0.97223 0.0625
        0.375 0.79904 0.86774
        0.3125 0.02777 0.9375
        0.0 0.35293 0.11774
        0.5 0.44446 0.375
        0.5 0.38975 0.10098
        0.125 0.79904 0.63227
        0.5 0.64707 0.61774
        0.625 0.20096 0.13227
        0.0 0.14707 0.61774
        0.1875 0.72306 0.66348
        0.0 0.88975 0.10098
        0.625 0.29904 0.63227
        0.6875 0.10373 0.5072
        0.8125 0.97223 0.4375
        0.6875 0.22306 0.66348
        0.6875 0.48528 0.8125
        0.1875 0.60373 0.5072
        0.75 0.15193 0.25
        0.1875 0.02777 0.5625
        0.5 0.46571 0.21872
        0.1875 0.32598 0.03122
        0.6875 0.02777 0.9375
        0.5 0.11026 0.60098
        0.625 0.09721 0.35098
        0.875 0.90279 0.85098
        0.8125 0.47223 0.0625
        0.8125 0.02777 0.5625
        0.8125 0.48528 0.6875
        0.1875 0.17402 0.53122
        0.0 0.64707 0.88227
        0.5 0.14707 0.88227
        0.0 0.81946 0.20196
        0.75 0.65193 0.25
        0.8125 0.67402 0.96879
        0.6875 0.47223 0.4375
        0.5 0.03429 0.71872
        0.625 0.90279 0.64902
        0.875 0.79904 0.63227
        0.125 0.40279 0.64902
        0.1875 0.47223 0.0625
        0.375 0.70096 0.36774
        0.875 0.20096 0.36774
        0.875 0.29904 0.86774
        0.1875 0.89627 0.0072
        0.6875 0.39627 0.0072
        0.8125 0.10373 0.99281
        0.3125 0.60373 0.99281
        0.375 0.59721 0.14902
        0.0 0.03429 0.78129
        0.3125 0.47223 0.4375
        0.6875 0.89627 0.49281
        0.1875 0.39627 0.49281
        0.25 0.65193 0.25
        0.125 0.09721 0.14902
        0.8125 0.89627 0.0072
        0.3125 0.39627 0.0072
        0.1875 0.10373 0.99281
        0.6875 0.60373 0.99281
        0.5 0.81946 0.29804
        0.6875 0.32598 0.46879
        0.375 0.40279 0.85098
        0.8125 0.98528 0.8125
        0.75 0.34808 0.75
        0.0 0.46571 0.28129
        0.3125 0.89627 0.49281
        0.8125 0.39627 0.49281
        0.5 0.18054 0.70196
        0.3125 0.98528 0.6875
        0.625 0.79904 0.86774
        0.6875 0.67402 0.53122
        0.125 0.29904 0.86774
        0.125 0.20096 0.36774
        0.3125 0.17402 0.96879
        0.5 0.61026 0.89902
        0.0 0.44446 0.125
        0.375 0.20096 0.13227
        0.375 0.29904 0.63227
        0.875 0.70096 0.13227
        0.0 0.53429 0.71872
        0.5 0.35293 0.38227
        0.3125 0.82598 0.03122
        0.8125 0.77695 0.16348
        0.8125 0.01473 0.1875
        0.3125 0.51473 0.1875
        0.25 0.34808 0.75
        0.1875 0.48528 0.6875
        0.5 0.85293 0.11774
        0.875 0.09721 0.14902
        0.0 0.94446 0.375
        0.8125 0.17402 0.53122
        0.5 0.53429 0.78129
        0.1875 0.67402 0.96879
        0.875 0.40279 0.64902
        0.6875 0.27695 0.16348
        0.1875 0.77695 0.16348
        0.5 0.68054 0.79804
        0.1875 0.01473 0.1875
        0.6875 0.51473 0.1875
        0.0 0.61026 0.60098
        0.125 0.59721 0.35098
        0.0 0.31946 0.29804
        0.6875 0.77695 0.33652
        0.5 0.88975 0.39902
        0.1875 0.82598 0.46879
        0.1875 0.27695 0.33652
        0.25 0.84808 0.75
        0.3125 0.72306 0.83652
        0.8125 0.22306 0.83652
        0.0 0.38975 0.39902
        0.3125 0.48528 0.8125
        0.5 0.96571 0.28129
        0.0 0.68054 0.70196
        0.3125 0.27695 0.16348
        """

        self.coord = "relative"

        self.cages = """
        14 0.5 0.91501 0.87486
        15 0.5 0.78114 0.59608
        15 0.0 0.21886 0.09608
        14 0.0 0.08499 0.37486
        12 0.0 0.0 0.0
        15 0.0 0.28114 0.59608
        12 -0.25 -0.61108 0.25
        14 0.0 0.58499 0.12514
        12 -0.25 -0.38892 -0.25
        14 0.5 0.41501 -0.37486
        14 0.5 1.08499 0.12514
        12 0.0 0.5 0.5
        12 -0.25 0.11108 0.75
        15 0.5 1.21886 0.40392
        12 0.25 0.38892 0.25
        14 0.0 0.41501 0.87486
        14 0.0 -0.08499 -0.37486
        12 -0.25 0.88892 0.25
        15 0.5 0.71886 0.09608
        14 0.5 0.58499 0.37486
        12 0.25 -0.88892 -0.25
        12 0.5 1.0 0.5
        15 0.5 0.28114 -0.09608
        12 0.25 -0.11108 -0.75
        15 0.0 -0.21886 -0.09608
        12 0.5 0.5 0.0
        15 0.0 0.71886 0.40392
        12 0.25 0.61108 -0.25
        """

        self.bondlen = 3

        self.cell = """
        13.518112406663548 36.32755830109184 17.45815101781426
        """

        self.density = 0.5578291804888631

        self.cell = cellvectors(a=13.518112406663548,
                                b=36.32755830109184,
                                c=17.45815101781426)
