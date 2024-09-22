# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (14,16,4,2,)
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
        9 70
        1 68
        14 163
        45 114
        99 66
        122 92
        52 87
        87 158
        167 74
        158 81
        117 50
        28 59
        117 111
        19 16
        34 40
        69 159
        143 23
        112 6
        166 64
        173 33
        124 108
        78 68
        112 100
        43 12
        132 149
        123 160
        192 156
        19 88
        192 82
        187 23
        152 196
        142 98
        38 119
        69 4
        202 15
        118 72
        46 92
        6 71
        181 3
        32 167
        5 181
        124 176
        113 191
        188 72
        141 59
        106 24
        27 168
        31 196
        83 204
        190 97
        177 147
        107 155
        62 134
        170 177
        168 23
        145 173
        200 71
        141 29
        180 105
        110 177
        2 48
        109 63
        49 97
        44 74
        109 70
        85 194
        10 96
        56 176
        55 5
        27 95
        84 173
        194 155
        138 182
        17 160
        202 135
        72 85
        178 61
        175 115
        99 53
        198 135
        165 48
        199 59
        201 98
        159 58
        66 205
        51 111
        186 86
        140 138
        164 174
        51 127
        65 123
        124 158
        29 200
        68 155
        109 178
        56 165
        119 114
        124 166
        20 67
        129 133
        157 16
        122 191
        35 94
        44 195
        2 67
        91 77
        8 74
        201 126
        69 88
        95 73
        61 151
        195 174
        104 18
        77 185
        79 48
        31 197
        57 5
        11 153
        166 84
        121 205
        121 48
        94 13
        101 46
        49 132
        115 80
        186 94
        44 131
        8 116
        169 60
        150 136
        42 120
        132 111
        21 18
        43 164
        20 195
        20 167
        24 162
        52 165
        202 111
        104 162
        139 66
        79 93
        25 50
        156 134
        54 29
        8 153
        92 183
        73 161
        81 200
        75 162
        59 158
        127 168
        189 172
        89 86
        179 62
        139 163
        155 15
        45 112
        22 91
        14 82
        147 76
        10 33
        73 145
        3 80
        130 60
        160 126
        84 39
        143 149
        159 46
        186 198
        21 58
        22 136
        71 80
        125 165
        17 98
        43 205
        25 76
        56 53
        19 35
        78 189
        179 83
        132 25
        22 37
        35 22
        31 182
        57 40
        11 75
        148 110
        90 120
        89 185
        101 4
        105 11
        55 103
        180 133
        131 192
        133 18
        41 104
        105 42
        127 149
        36 15
        19 182
        154 184
        50 86
        23 161
        181 110
        41 184
        200 128
        188 196
        73 10
        161 3
        194 135
        1 188
        90 172
        159 191
        138 184
        63 205
        180 65
        131 126
        100 204
        127 86
        83 141
        129 126
        145 103
        146 106
        143 57
        113 154
        89 13
        102 95
        96 147
        2 62
        24 123
        44 123
        203 17
        0 156
        30 183
        56 204
        144 70
        52 171
        197 15
        149 185
        43 7
        102 143
        83 137
        97 76
        103 3
        77 187
        21 42
        102 25
        79 63
        131 47
        177 145
        178 66
        37 51
        189 58
        118 90
        32 163
        40 115
        32 63
        5 185
        91 168
        7 192
        28 79
        199 137
        36 186
        99 199
        191 172
        109 67
        27 50
        93 163
        128 115
        128 114
        87 54
        202 136
        96 26
        1 88
        36 117
        150 197
        178 137
        61 156
        154 85
        100 9
        27 190
        193 184
        180 201
        121 53
        38 87
        105 142
        125 100
        153 17
        203 172
        139 12
        31 68
        117 152
        38 60
        2 144
        38 64
        78 4
        106 74
        157 198
        154 69
        108 204
        67 151
        129 92
        116 98
        142 183
        81 84
        107 72
        175 26
        179 53
        189 30
        199 176
        12 82
        108 29
        39 71
        26 57
        101 30
        125 54
        187 181
        89 76
        108 6
        170 95
        116 47
        122 160
        41 101
        171 64
        182 4
        55 147
        60 175
        41 120
        122 162
        85 157
        94 197
        141 144
        1 58
        120 75
        14 174
        81 119
        161 40
        148 169
        36 140
        193 90
        49 77
        203 75
        125 144
        97 170
        148 80
        99 93
        107 193
        152 150
        7 151
        0 47
        49 136
        113 104
        201 82
        65 174
        62 93
        78 193
        138 157
        170 187
        110 26
        42 30
        12 116
        146 142
        130 128
        175 33
        35 135
        130 54
        179 61
        91 13
        139 0
        190 13
        34 10
        46 18
        55 34
        176 171
        146 129
        140 196
        8 164
        65 153
        20 134
        0 167
        119 33
        6 169
        37 152
        188 16
        9 171
        9 137
        133 24
        190 150
        14 134
        45 39
        28 70
        51 198
        169 166
        118 113
        194 88
        118 21
        106 11
        146 47
        121 151
        34 114
        103 39
        148 173
        140 107
        32 164
        7 195
        28 52
        37 16
        203 183
        112 130
        102 96
        45 64
        """

        self.waters = """
        0.875 0.0 0.83736
        0.81645 0.5 0.13573
        0.18355 0.30855 0.73106
        0.81772 0.18591 0.45772
        0.75 0.0 0.11413
        0.5 0.18591 0.41648
        0.69272 0.31091 0.56923
        0.5 0.375 0.8251
        0.375 0.81772 0.88722
        0.5 0.0 0.66098
        0.18228 0.81409 0.45772
        0.18591 0.81772 0.95772
        0.69145 0.69272 0.84901
        0.69272 0.0 0.30481
        0.0 0.5 0.83902
        0.5 0.69145 0.22076
        0.0 0.30728 0.19519
        0.5 0.625 0.95907
        0.0 0.30728 0.02829
        0.81645 0.18228 0.18655
        0.18355 0.18228 0.81345
        0.0 0.5 0.05658
        0.81645 0.30855 0.26894
        0.0 0.25 0.38587
        0.18591 0.18228 0.95772
        0.30728 0.69145 0.34901
        0.375 0.5 0.45907
        0.0 0.875 0.33736
        0.18228 0.81645 0.68655
        0.0 0.25 0.61413
        0.81409 0.81772 0.04228
        0.81645 0.81772 0.18655
        0.18355 0.81772 0.81345
        0.125 0.68909 0.5
        0.30728 0.0 0.47171
        0.69145 0.18355 0.23106
        0.30855 0.81645 0.23106
        0.0 0.375 0.24834
        0.30728 0.68909 0.56923
        0.69272 0.0 0.52829
        0.18228 0.18591 0.45772
        0.125 0.0 0.07057
        0.0 0.69272 0.02829
        0.5 0.625 0.8251
        0.375 0.18228 0.88722
        0.5 0.0 0.55658
        0.81409 0.18228 0.04228
        0.75 0.0 0.88587
        0.30855 0.5 0.72076
        0.625 0.5 0.3251
        0.18228 0.81645 0.31345
        0.18355 0.30855 0.26894
        0.30728 0.69145 0.65099
        0.69145 0.5 0.72076
        0.18228 0.375 0.61278
        0.5 0.0 0.44342
        0.625 0.5 0.6749
        0.30728 0.31091 0.43078
        0.81409 0.5 0.08353
        0.0 0.875 0.66264
        0.375 0.5 0.54094
        0.69145 0.18355 0.76894
        0.0 0.375 0.75166
        0.30855 0.81645 0.76894
        0.5 0.81409 0.58353
        0.18591 0.5 0.91648
        0.69145 0.81645 0.76894
        0.30855 0.18355 0.76894
        0.69145 0.69272 0.15099
        0.625 0.18228 0.11278
        0.30728 0.0 0.69519
        0.81772 0.18591 0.54228
        0.18355 0.5 0.13573
        0.0 0.875 0.42943
        0.25 0.0 0.88587
        0.31091 0.875 0.0
        0.5 0.81645 0.36428
        0.69272 0.30855 0.34901
        0.625 0.81772 0.11278
        0.18355 0.69145 0.73106
        0.875 0.31091 0.5
        0.0 0.875 0.57057
        0.81645 0.5 0.86428
        0.81772 0.18355 0.68655
        0.81772 0.81409 0.54228
        0.30855 0.30728 0.15099
        0.30728 0.0 0.30481
        0.18228 0.625 0.61278
        0.69145 0.30728 0.15099
        0.5 0.0 0.33902
        0.31091 0.69272 0.06923
        0.81772 0.18355 0.31345
        0.68909 0.125 0.0
        0.0 0.625 0.75166
        0.625 0.0 0.25166
        0.0 0.75 0.38587
        0.30728 0.68909 0.43078
        0.69272 0.69145 0.34901
        0.68909 0.69272 0.93078
        0.81645 0.69145 0.73106
        0.5 0.18355 0.63573
        0.875 0.0 0.07057
        0.18228 0.625 0.38722
        0.69272 0.0 0.47171
        0.18591 0.18228 0.04228
        0.0 0.69272 0.97171
        0.125 0.0 0.92943
        0.30855 0.69272 0.15099
        0.81772 0.375 0.61278
        0.375 0.0 0.74834
        0.625 0.5 0.45907
        0.30855 0.5 0.27924
        0.5 0.18591 0.58353
        0.31091 0.30728 0.06923
        0.30728 0.0 0.52829
        0.125 0.31091 0.5
        0.625 0.81772 0.88722
        0.18355 0.69145 0.26894
        0.18591 0.5 0.08353
        0.18228 0.81409 0.54228
        0.18591 0.81772 0.04228
        0.5 0.5 0.75
        0.5 0.25 0.0
        0.31091 0.30728 0.93078
        0.81772 0.625 0.61278
        0.30728 0.30855 0.65099
        0.68909 0.30728 0.93078
        0.18228 0.18355 0.31345
        0.18228 0.18591 0.54228
        0.81409 0.18228 0.95772
        0.30728 0.31091 0.56923
        0.625 0.18228 0.88722
        0.375 0.5 0.3251
        0.0 0.30728 0.97171
        0.0 0.30728 0.80481
        0.5 0.30855 0.22076
        0.69145 0.5 0.27924
        0.69272 0.0 0.69519
        0.125 0.0 0.16264
        0.81645 0.81772 0.81345
        0.18355 0.81772 0.18655
        0.0 0.125 0.66264
        0.81409 0.81772 0.95772
        0.18228 0.375 0.38722
        0.18228 0.18355 0.68655
        0.81772 0.81409 0.45772
        0.875 0.0 0.92943
        0.5 0.81409 0.41648
        0.75 0.5 0.5
        0.30728 0.30855 0.34901
        0.81645 0.69145 0.26894
        0.5 0.30855 0.77924
        0.0 0.625 0.24834
        0.31091 0.69272 0.93078
        0.375 0.18228 0.11278
        0.5 0.625 0.1749
        0.81645 0.18228 0.81345
        0.18355 0.18228 0.18655
        0.0 0.75 0.61413
        0.68909 0.30728 0.06923
        0.5 0.375 0.95907
        0.0 0.125 0.42943
        0.31091 0.125 0.0
        0.0 0.69272 0.80481
        0.30855 0.69272 0.84901
        0.375 0.5 0.6749
        0.69272 0.68909 0.56923
        0.125 0.0 0.83736
        0.0 0.125 0.33736
        0.625 0.5 0.54094
        0.81772 0.625 0.38722
        0.5 0.81645 0.63573
        0.5 0.625 0.04094
        0.875 0.68909 0.5
        0.18355 0.5 0.86428
        0.25 0.5 0.5
        0.69272 0.69145 0.65099
        0.69272 0.68909 0.43078
        0.625 0.0 0.74834
        0.81645 0.30855 0.73106
        0.0 0.5 0.94342
        0.69272 0.31091 0.43078
        0.875 0.0 0.16264
        0.68909 0.875 0.0
        0.25 0.0 0.11413
        0.5 0.18355 0.36428
        0.375 0.0 0.25166
        0.81772 0.375 0.38722
        0.0 0.5 0.16098
        0.68909 0.69272 0.06923
        0.81772 0.81645 0.31345
        0.5 0.375 0.04094
        0.69145 0.30728 0.84901
        0.375 0.81772 0.11278
        0.5 0.375 0.1749
        0.30855 0.30728 0.84901
        0.0 0.69272 0.19519
        0.69145 0.81645 0.23106
        0.30855 0.18355 0.23106
        0.81772 0.81645 0.68655
        0.0 0.125 0.57057
        0.81409 0.5 0.91648
        0.5 0.5 0.25
        0.5 0.75 0.0
        0.69272 0.30855 0.65099
        0.5 0.69145 0.77924
        """

        self.coord = "relative"

        self.cages = """
        12 0.0 0.0 0.5
        14 0.0 0.5 0.55596
        12 -0.2342 0.5 -0.2088
        15 0.5 0.0 0.17423
        14 0.0 0.5 -0.55596
        15 -0.5 0.0 -0.17423
        14 -0.5 0.24364 0.5
        12 0.2342 0.5 0.2088
        14 0.0 -0.22913 0.11316
        14 -0.22913 0.0 0.61316
        12 0.0 0.0 0.0
        12 -0.5 0.2342 0.7088
        14 0.0 -0.22913 -0.11316
        12 0.5 0.5 -0.60778
        14 -0.22913 0.0 0.38684
        12 0.2342 -0.5 -0.2088
        15 0.0 0.5 -0.67423
        15 0.0 0.5 0.67423
        14 0.0 0.22913 -0.11316
        14 0.5 -0.24364 0.5
        16 0.0 0.0 0.75
        14 0.22913 0.0 0.61316
        14 0.0 0.22913 0.11316
        14 0.5 0.0 -1.05596
        12 -0.5 0.5 -0.10778
        12 -0.2342 -0.5 0.2088
        16 0.0 0.0 0.25
        14 0.24364 0.5 0.0
        12 0.5 0.5 0.10778
        14 -0.5 0.0 1.05596
        14 -0.24364 -0.5 0.0
        12 0.5 0.2342 -0.7088
        12 0.5 -0.2342 -0.7088
        14 0.22913 0.0 -0.61316
        12 -0.5 0.5 0.60778
        12 0.5 -0.2342 0.7088
        """

        self.bondlen = 3

        self.cell = """
        13.002551975998497 13.002551975998497 57.158147048395534
        """

        self.density = 0.6371822012607398

        self.cell = cellvectors(a=13.002551975998497,
                                b=13.002551975998497,
                                c=57.158147048395534)