# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (10,16,4,0,)
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
        91 5
        9 71
        99 58
        62 7
        118 78
        136 100
        135 103
        2 73
        102 68
        146 126
        83 133
        119 93
        106 82
        33 78
        33 102
        125 90
        23 6
        62 87
        155 27
        79 95
        71 63
        24 81
        143 159
        151 32
        115 66
        1 10
        0 120
        170 77
        116 114
        90 124
        118 157
        79 114
        167 43
        18 67
        15 123
        114 67
        150 27
        18 75
        92 69
        125 20
        121 141
        169 26
        34 38
        23 103
        149 142
        51 80
        111 3
        107 40
        132 6
        65 73
        169 139
        170 141
        137 160
        86 58
        33 155
        60 30
        8 67
        26 54
        43 67
        41 131
        119 11
        52 51
        35 66
        112 163
        113 92
        47 158
        4 143
        102 70
        170 74
        13 55
        103 85
        0 76
        136 81
        95 135
        48 124
        63 92
        37 27
        9 45
        168 64
        8 49
        127 7
        94 111
        125 155
        48 28
        15 110
        61 69
        130 97
        89 104
        11 160
        128 138
        54 82
        105 32
        83 7
        86 120
        158 139
        148 144
        36 164
        163 49
        132 49
        94 88
        13 113
        123 32
        117 153
        171 60
        146 52
        57 148
        73 143
        57 121
        128 31
        14 12
        145 5
        21 165
        14 171
        149 50
        54 98
        68 109
        147 80
        76 109
        17 162
        148 129
        115 151
        162 16
        38 138
        82 165
        34 168
        70 5
        163 85
        139 69
        8 39
        22 46
        146 71
        4 46
        144 88
        127 10
        101 162
        46 150
        29 44
        91 86
        36 137
        153 29
        0 60
        22 94
        128 24
        30 109
        100 165
        112 116
        140 49
        53 115
        2 122
        107 96
        135 75
        65 74
        96 56
        55 167
        61 138
        94 143
        107 44
        72 166
        126 106
        47 13
        90 17
        127 72
        160 42
        166 121
        36 132
        35 130
        39 167
        10 73
        20 145
        59 43
        2 134
        64 45
        41 140
        137 110
        116 113
        118 90
        6 32
        52 31
        1 77
        152 169
        37 17
        41 84
        118 145
        34 31
        47 98
        125 16
        63 59
        155 29
        99 124
        12 76
        122 74
        134 166
        55 38
        135 26
        156 70
        53 84
        25 19
        59 98
        64 100
        133 101
        110 151
        156 86
        21 126
        62 111
        136 171
        14 5
        4 7
        2 50
        124 29
        33 89
        142 19
        119 23
        113 43
        91 60
        9 51
        57 133
        79 164
        20 102
        72 122
        36 161
        107 150
        98 75
        129 154
        79 108
        13 64
        153 111
        152 108
        89 44
        157 104
        80 30
        92 45
        112 69
        22 148
        53 105
        137 19
        37 96
        44 58
        147 24
        8 95
        57 87
        4 101
        169 85
        65 88
        147 120
        95 6
        93 108
        82 139
        38 165
        168 51
        157 99
        59 106
        158 100
        16 58
        131 142
        61 106
        40 87
        103 140
        37 78
        56 99
        119 164
        156 80
        20 120
        153 40
        52 70
        104 30
        117 22
        160 84
        1 19
        25 123
        168 171
        130 50
        42 93
        127 170
        34 12
        129 101
        42 18
        156 89
        157 91
        61 167
        23 131
        164 151
        161 114
        25 132
        116 152
        117 154
        31 68
        63 158
        110 77
        88 141
        50 105
        56 28
        149 66
        117 27
        11 142
        14 146
        10 149
        130 144
        150 16
        87 46
        108 75
        41 25
        3 83
        83 166
        35 121
        159 134
        26 39
        28 40
        97 123
        131 105
        55 54
        48 3
        65 97
        136 71
        21 81
        112 39
        21 12
        72 66
        47 152
        154 48
        28 133
        159 3
        128 126
        154 17
        35 77
        1 97
        78 109
        129 159
        161 163
        45 138
        15 74
        0 81
        9 24
        56 162
        53 122
        93 85
        134 144
        62 141
        145 76
        11 115
        15 84
        96 104
        161 42
        147 68
        140 18
        """

        self.waters="""
        0.5 0.73043 0.80821
        0.5 0.26957 0.19179
        0.8125 0.31867 0.80821
        0.0 0.45093 0.875
        0.5 0.43296 0.0
        0.0 0.71031 0.875
        0.75 0.13699 0.5
        0.31369 0.40074 0.0
        0.81369 0.05038 0.31679
        0.3125 0.81867 0.19179
        0.5 0.30796 0.0
        0.3125 0.18133 0.80821
        0.8125 0.77898 0.68321
        0.125 0.9161 0.5
        0.0 0.76957 0.80821
        0.0 0.23043 0.19179
        0.5 0.56705 0.0
        0.68632 0.55038 0.68321
        0.0 0.06705 0.0
        0.5 0.21031 0.125
        0.6875 0.66176 0.0
        0.6875 0.81867 0.80821
        0.68632 0.44963 0.31679
        0.625 0.13569 0.68321
        0.5 0.78969 0.125
        0.6875 0.18133 0.19179
        0.68869 0.0 0.625
        0.68632 0.55038 0.31679
        0.18869 0.5 0.625
        0.0 0.54907 0.125
        0.375 0.69205 0.5
        0.8125 0.77898 0.31679
        0.875 0.19205 0.5
        0.875 0.63569 0.31679
        0.875 0.80796 0.5
        0.375 0.30796 0.5
        0.375 0.13569 0.31679
        0.625 0.5839 0.5
        0.75 0.86301 0.5
        0.68869 0.0 0.375
        0.18869 0.5 0.375
        0.8125 0.16176 0.0
        0.18632 0.09926 0.0
        0.0 0.96648 0.19179
        0.18869 0.58259 0.19179
        0.375 0.86432 0.31679
        0.5 0.46648 0.19179
        0.18632 0.94963 0.68321
        0.0 0.5 0.75
        0.68869 0.08259 0.19179
        0.6875 0.27898 0.68321
        0.1875 0.77898 0.31679
        0.0 0.76957 0.19179
        0.0 0.23043 0.80821
        0.81369 0.94963 0.68321
        0.875 0.9161 0.5
        0.31369 0.55038 0.68321
        0.375 0.4161 0.5
        0.31369 0.59926 0.0
        0.0 0.93296 0.0
        0.3125 0.72102 0.68321
        0.68869 0.91741 0.19179
        0.18869 0.41741 0.19179
        0.18632 0.90074 0.0
        0.25 0.86301 0.5
        0.8125 0.31867 0.19179
        0.3125 0.27898 0.68321
        0.0 0.03352 0.19179
        0.6875 0.72102 0.31679
        0.5 0.95093 0.125
        0.0 0.71031 0.125
        0.1875 0.83824 0.0
        0.1875 0.31867 0.80821
        0.6875 0.33824 0.0
        0.0 0.28969 0.125
        0.0 0.03352 0.80821
        0.6875 0.72102 0.68321
        0.3125 0.27898 0.31679
        0.75 0.63699 0.5
        0.125 0.0839 0.5
        0.3125 0.72102 0.31679
        0.5 0.78969 0.875
        0.68869 0.91741 0.80821
        0.18869 0.41741 0.80821
        0.0 0.19205 0.0
        0.5 0.04907 0.875
        0.3125 0.66176 0.0
        0.31369 0.44963 0.31679
        0.875 0.36432 0.31679
        0.125 0.63569 0.31679
        0.81132 0.58259 0.80821
        0.1875 0.68133 0.80821
        0.31132 0.91741 0.19179
        0.31132 0.08259 0.80821
        0.81132 0.41741 0.19179
        0.875 0.0839 0.5
        0.375 0.5839 0.5
        0.6875 0.27898 0.31679
        0.0 0.96648 0.80821
        0.18869 0.58259 0.80821
        0.375 0.86432 0.68321
        0.5 0.46648 0.80821
        0.8125 0.68133 0.19179
        0.68869 0.08259 0.80821
        0.25 0.63699 0.5
        0.8125 0.22102 0.68321
        0.81369 0.90074 0.0
        0.31369 0.55038 0.31679
        0.18632 0.05038 0.68321
        0.625 0.69205 0.5
        0.1875 0.22102 0.31679
        0.0 0.45093 0.125
        0.5 0.0 0.25
        0.18632 0.94963 0.31679
        0.18632 0.05038 0.31679
        0.1875 0.22102 0.68321
        0.31132 0.0 0.375
        0.81132 0.5 0.375
        0.875 0.63569 0.68321
        0.375 0.13569 0.68321
        0.5 0.69205 0.0
        0.25 0.36301 0.5
        0.0 0.28969 0.875
        0.8125 0.22102 0.31679
        0.0 0.54907 0.875
        0.68632 0.59926 0.0
        0.8125 0.83824 0.0
        0.3125 0.33824 0.0
        0.6875 0.81867 0.19179
        0.68632 0.44963 0.68321
        0.625 0.30796 0.5
        0.6875 0.18133 0.80821
        0.625 0.13569 0.31679
        0.31369 0.44963 0.68321
        0.875 0.36432 0.68321
        0.81369 0.05038 0.68321
        0.3125 0.81867 0.80821
        0.3125 0.18133 0.19179
        0.625 0.86432 0.31679
        0.5 0.95093 0.875
        0.81369 0.09926 0.0
        0.125 0.36432 0.31679
        0.5 0.21031 0.875
        0.68632 0.40074 0.0
        0.75 0.36301 0.5
        0.8125 0.68133 0.80821
        0.0 0.80796 0.0
        0.5 0.73043 0.19179
        0.625 0.4161 0.5
        0.5 0.26957 0.80821
        0.5 0.53352 0.19179
        0.125 0.19205 0.5
        0.31132 0.0 0.625
        0.0 0.5 0.25
        0.81132 0.5 0.625
        0.81132 0.58259 0.19179
        0.1875 0.68133 0.19179
        0.125 0.63569 0.68321
        0.31132 0.91741 0.80821
        0.81132 0.41741 0.80821
        0.1875 0.16176 0.0
        0.31132 0.08259 0.19179
        0.5 0.53352 0.80821
        0.5 0.04907 0.125
        0.25 0.13699 0.5
        0.625 0.86432 0.68321
        0.125 0.36432 0.68321
        0.81369 0.94963 0.31679
        0.125 0.80796 0.5
        0.5 0.0 0.75
        0.1875 0.31867 0.19179
        0.1875 0.77898 0.68321
        """

        self.coord= "relative"

        self.cages="""
        12 0.0 0.0 0.5
        14 0.5 0.63409 0.26717
        14 1.0 0.56741 0.5
        15 0.5 0.78762 0.5
        15 0.0 -0.28762 0.5
        14 0.0 -0.13409 0.26717
        14 -0.5 -0.06741 0.5
        12 0.5 0.5 0.5
        12 0.25 -0.25 0.0
        14 0.5 0.36591 -0.26717
        14 -0.24526 0.0 0.0
        12 0.5 0.87114 0.0
        14 0.0 0.13409 -0.26717
        14 0.5 0.63409 -0.26717
        12 0.25 0.25 0.0
        14 0.5 0.06741 0.5
        14 0.74526 0.5 0.0
        15 0.0 0.28762 0.5
        12 0.5 0.12886 0.0
        14 0.25474 0.5 0.0
        14 0.5 0.36591 0.26717
        14 0.0 0.43259 0.5
        12 -0.25 0.25 0.0
        14 0.0 -0.13409 -0.26717
        14 0.24526 0.0 0.0
        12 0.0 0.37114 0.0
        12 -0.25 -0.25 0.0
        12 0.0 -0.37114 0.0
        15 0.5 0.21238 0.5
        14 0.0 0.13409 0.26717
        """

        self.bondlen = 3


        self.cell = """
        12.85918946354857 48.068047982839694 12.990393566570004
        """

        self.density = 0.6402768662754011



        self.cell = cellvectors(a=12.85918946354857,
                           b=48.068047982839694,
                           c=12.990393566570004)
