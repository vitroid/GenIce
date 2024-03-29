# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (20,8,8,4,)
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
        100 33
        127 89
        161 16
        183 144
        148 107
        37 28
        91 6
        173 17
        34 53
        191 69
        16 216
        154 71
        156 72
        183 55
        44 226
        81 177
        198 137
        198 136
        42 168
        74 195
        143 132
        35 25
        219 12
        146 98
        158 177
        224 30
        114 156
        154 113
        96 179
        180 140
        83 153
        70 92
        6 196
        110 65
        55 2
        7 208
        194 206
        142 108
        129 227
        90 16
        109 82
        49 177
        41 40
        81 208
        112 135
        169 138
        97 17
        78 53
        197 19
        88 118
        64 72
        52 223
        179 24
        180 26
        54 59
        121 173
        58 176
        132 199
        159 10
        15 43
        93 62
        164 151
        221 165
        11 209
        114 87
        85 122
        85 123
        139 186
        188 57
        63 4
        138 31
        138 32
        181 166
        13 97
        1 177
        7 181
        110 74
        94 23
        45 59
        125 66
        73 136
        115 190
        28 193
        83 25
        119 117
        9 11
        82 98
        133 157
        10 6
        202 107
        203 107
        190 105
        20 161
        36 18
        21 162
        50 29
        125 86
        204 41
        212 59
        149 117
        105 72
        65 214
        204 60
        3 40
        185 3
        39 124
        146 81
        8 192
        5 188
        126 201
        111 67
        167 130
        141 213
        27 149
        120 209
        121 211
        79 157
        83 207
        48 19
        111 202
        40 2
        139 62
        42 166
        145 119
        99 18
        216 205
        198 33
        0 12
        54 216
        66 227
        4 160
        64 225
        38 118
        178 105
        225 210
        31 165
        33 166
        182 94
        167 105
        155 209
        109 106
        112 104
        49 98
        11 131
        32 17
        101 193
        46 70
        121 223
        81 205
        34 108
        222 2
        158 189
        15 68
        22 224
        194 210
        224 76
        115 143
        37 156
        120 89
        134 57
        50 102
        82 54
        97 211
        43 31
        202 190
        129 34
        19 187
        175 61
        192 217
        226 18
        196 217
        52 89
        115 203
        148 67
        121 57
        55 75
        159 161
        50 170
        52 172
        27 73
        169 221
        118 106
        96 9
        141 44
        46 136
        39 157
        113 217
        150 134
        160 61
        114 218
        146 59
        147 60
        133 193
        147 119
        116 29
        172 206
        220 35
        148 137
        221 36
        197 144
        126 212
        158 80
        15 160
        127 213
        224 182
        51 122
        78 220
        20 46
        162 187
        52 134
        46 77
        9 150
        152 5
        227 71
        111 79
        156 130
        169 173
        58 35
        102 92
        14 71
        146 122
        32 182
        3 75
        70 137
        69 135
        106 49
        9 56
        47 131
        100 3
        185 145
        165 76
        77 166
        85 6
        116 205
        73 100
        139 124
        99 206
        67 197
        95 24
        141 93
        92 187
        153 88
        38 207
        104 60
        16 151
        63 68
        47 176
        37 191
        184 122
        67 143
        203 70
        180 40
        62 135
        95 32
        137 19
        184 45
        219 193
        38 200
        125 97
        128 79
        36 186
        141 169
        91 54
        176 211
        110 210
        55 29
        144 222
        0 101
        152 155
        161 170
        120 5
        88 24
        69 72
        170 212
        172 213
        117 2
        84 212
        203 199
        66 217
        60 107
        8 1
        184 100
        50 201
        130 222
        76 96
        189 14
        225 214
        140 77
        151 14
        45 75
        86 113
        188 160
        74 43
        214 69
        22 221
        222 87
        126 75
        138 30
        148 26
        85 192
        20 102
        140 117
        127 173
        149 29
        189 207
        163 209
        164 208
        39 61
        142 182
        13 35
        219 99
        96 200
        199 140
        41 145
        27 7
        36 93
        57 56
        25 24
        62 65
        189 106
        191 104
        63 210
        215 114
        79 101
        12 171
        219 191
        133 112
        109 196
        88 155
        198 168
        213 61
        168 48
        185 201
        47 163
        211 103
        223 94
        133 175
        130 119
        14 196
        109 153
        39 194
        87 41
        104 171
        17 23
        13 53
        183 48
        115 174
        142 200
        0 64
        108 118
        164 123
        128 218
        225 186
        44 65
        135 218
        20 42
        4 206
        179 31
        33 180
        7 116
        53 66
        51 73
        128 64
        159 84
        167 174
        68 186
        116 45
        8 164
        175 99
        162 145
        195 76
        179 152
        28 215
        47 227
        28 132
        82 80
        158 216
        83 80
        89 23
        90 91
        147 199
        220 80
        21 197
        51 10
        153 113
        202 171
        110 124
        123 90
        22 226
        185 149
        163 23
        174 92
        1 129
        175 93
        226 172
        176 94
        1 220
        204 26
        112 215
        159 136
        44 30
        178 218
        25 200
        78 98
        195 150
        43 5
        95 108
        152 56
        131 103
        13 95
        142 58
        86 103
        127 188
        68 165
        194 12
        139 128
        84 91
        111 178
        22 223
        37 147
        18 214
        204 215
        126 48
        181 123
        15 124
        74 30
        27 77
        129 58
        38 11
        170 205
        134 4
        42 90
        201 187
        101 143
        63 195
        157 171
        183 102
        78 192
        132 26
        34 49
        103 56
        84 168
        125 163
        181 184
        120 150
        21 87
        10 151
        0 190
        174 144
        131 154
        51 208
        207 154
        8 71
        167 162
        21 178
        155 86
        """

        self.waters = """
        0.19506 0.78232 0.9375
        0.69506 0.28232 0.9375
        0.0 0.37223 0.75
        0.94805 0.30836 0.4375
        0.35757 0.89766 0.875
        0.44805 0.80836 0.4375
        0.77043 0.91346 0.4375
        0.85416 0.24034 0.9375
        0.71549 0.11887 0.875
        0.5 0.77679 0.9375
        0.80154 0.95191 0.25
        0.55195 0.80836 0.0625
        0.22957 0.91346 0.0625
        0.58868 0.31664 0.5625
        0.71549 0.88114 0.125
        0.35416 0.74034 0.5625
        0.80495 0.78232 0.9375
        0.5 0.22321 0.4375
        0.30495 0.28232 0.9375
        0.0 0.72321 0.4375
        0.91132 0.81664 0.9375
        0.08527 0.56389 0.4375
        0.41132 0.31664 0.9375
        0.5 0.12777 0.25
        0.55535 0.56112 0.5625
        0.59208 0.53155 0.75
        0.05535 0.06112 0.5625
        0.91132 0.18336 0.0625
        0.14924 0.14766 0.75
        0.91473 0.43612 0.9375
        0.40792 0.46845 0.25
        0.44465 0.56112 0.5625
        0.5 0.37777 0.375
        0.94465 0.06112 0.5625
        0.64924 0.35234 0.25
        0.60679 0.37776 0.75
        0.32538 0.36887 0.75
        0.14584 0.24034 0.9375
        0.58868 0.68336 0.0625
        0.28452 0.88114 0.375
        0.0 0.27679 0.5625
        0.05195 0.30836 0.4375
        0.89322 0.87776 0.75
        0.41132 0.68336 0.4375
        0.35076 0.35234 0.25
        0.85757 0.39766 0.625
        0.94465 0.93889 0.0625
        0.58527 0.06389 0.0625
        0.94805 0.69164 0.5625
        0.69846 0.45191 0.25
        0.91473 0.56389 0.0625
        0.82538 0.13114 0.25
        0.41473 0.06389 0.0625
        0.64584 0.25966 0.4375
        0.78452 0.61887 0.625
        0.94805 0.46844 0.75
        0.5 0.87223 0.75
        0.44805 0.96844 0.75
        0.58868 0.31664 0.9375
        0.81903 0.5 0.5
        0.10679 0.12224 0.25
        0.31903 0.0 0.5
        0.27043 0.41346 0.4375
        0.35416 0.74034 0.9375
        0.21549 0.61887 0.875
        0.30154 0.45191 0.25
        0.64243 0.10234 0.375
        0.08868 0.81664 0.5625
        0.35076 0.64766 0.75
        0.21549 0.38114 0.125
        0.0 0.87777 0.125
        0.68097 0.0 0.0
        0.18097 0.5 0.0
        0.89322 0.12224 0.25
        0.39322 0.62224 0.25
        0.91473 0.43612 0.5625
        0.44465 0.56112 0.9375
        0.94465 0.06112 0.9375
        0.69506 0.28232 0.5625
        0.19506 0.78232 0.5625
        0.69846 0.54809 0.75
        0.78452 0.38114 0.125
        0.72957 0.58655 0.5625
        0.64924 0.64766 0.75
        0.85416 0.75966 0.4375
        0.77043 0.08655 0.5625
        0.58527 0.93612 0.5625
        0.08527 0.43612 0.5625
        0.58868 0.68336 0.4375
        0.44805 0.03157 0.25
        0.82538 0.86887 0.75
        0.80495 0.78232 0.5625
        0.0 0.72321 0.0625
        0.30495 0.28232 0.5625
        0.5 0.22321 0.0625
        0.55535 0.43889 0.4375
        0.5 0.62223 0.875
        0.55195 0.19164 0.5625
        0.72957 0.41346 0.4375
        0.28452 0.11887 0.875
        0.91132 0.18336 0.4375
        0.17463 0.86887 0.75
        0.94805 0.69164 0.9375
        0.55195 0.96844 0.75
        0.17463 0.13114 0.25
        0.14243 0.60234 0.125
        0.67463 0.63114 0.25
        0.09208 0.96845 0.25
        0.59208 0.46845 0.25
        0.69506 0.71768 0.4375
        0.32538 0.63114 0.25
        0.14584 0.75966 0.4375
        0.19506 0.21768 0.4375
        0.64243 0.89766 0.625
        0.14243 0.39766 0.625
        0.08868 0.81664 0.9375
        0.85757 0.39766 0.875
        0.0 0.27679 0.9375
        0.60679 0.62224 0.25
        0.05195 0.30836 0.0625
        0.46668 0.87224 0.25
        0.46668 0.12776 0.75
        0.80495 0.21768 0.4375
        0.80154 0.04809 0.75
        0.30495 0.71768 0.4375
        0.58527 0.06389 0.4375
        0.91473 0.56389 0.4375
        0.41473 0.06389 0.4375
        0.21549 0.61887 0.625
        0.64584 0.25966 0.0625
        0.08527 0.43612 0.9375
        0.58527 0.93612 0.9375
        0.09208 0.03155 0.75
        0.22957 0.08655 0.5625
        0.41473 0.93612 0.9375
        0.21549 0.38114 0.375
        0.90792 0.96845 0.25
        0.0 0.87777 0.375
        0.44465 0.43889 0.4375
        0.27043 0.58655 0.5625
        0.0 0.12223 0.875
        0.35416 0.25966 0.4375
        0.55535 0.43889 0.0625
        0.10679 0.87776 0.75
        0.03333 0.62776 0.75
        0.03333 0.37224 0.25
        0.78452 0.38114 0.375
        0.08868 0.18336 0.0625
        0.05535 0.93889 0.4375
        0.94805 0.30836 0.0625
        0.44805 0.80836 0.0625
        0.77043 0.91346 0.0625
        0.5 0.77679 0.5625
        0.64584 0.74034 0.5625
        0.64243 0.89766 0.875
        0.55195 0.80836 0.4375
        0.14243 0.39766 0.875
        0.22957 0.91346 0.4375
        0.72957 0.58655 0.9375
        0.85076 0.85234 0.25
        0.35757 0.89766 0.625
        0.85416 0.75966 0.0625
        0.05195 0.53157 0.25
        0.55195 0.03157 0.25
        0.77043 0.08655 0.9375
        0.40792 0.53155 0.75
        0.90792 0.03155 0.75
        0.08527 0.56389 0.0625
        0.91132 0.81664 0.5625
        0.41132 0.31664 0.5625
        0.85757 0.60234 0.125
        0.19846 0.95191 0.25
        0.35757 0.10234 0.125
        0.44805 0.19164 0.5625
        0.05195 0.69164 0.9375
        0.28452 0.11887 0.625
        0.55195 0.19164 0.9375
        0.72957 0.41346 0.0625
        0.14243 0.60234 0.375
        0.5 0.62223 0.625
        0.0 0.12223 0.625
        0.85076 0.14766 0.75
        0.5 0.37777 0.125
        0.96668 0.62776 0.75
        0.85416 0.24034 0.5625
        0.96668 0.37224 0.25
        0.30154 0.54809 0.75
        0.0 0.62777 0.25
        0.41473 0.93612 0.5625
        0.69506 0.71768 0.0625
        0.14584 0.75966 0.0625
        0.19506 0.21768 0.0625
        0.71549 0.11887 0.625
        0.19846 0.04809 0.75
        0.28452 0.88114 0.125
        0.41132 0.68336 0.0625
        0.71549 0.88114 0.375
        0.05195 0.69164 0.5625
        0.94465 0.93889 0.4375
        0.05535 0.06112 0.9375
        0.55535 0.56112 0.9375
        0.94805 0.53157 0.25
        0.14924 0.85234 0.25
        0.05535 0.93889 0.0625
        0.08868 0.18336 0.4375
        0.81903 0.5 0.0
        0.31903 0.0 0.0
        0.64584 0.74034 0.9375
        0.80495 0.21768 0.0625
        0.53333 0.87224 0.25
        0.30495 0.71768 0.0625
        0.53333 0.12776 0.75
        0.85757 0.60234 0.375
        0.35757 0.10234 0.375
        0.27043 0.41346 0.0625
        0.14584 0.24034 0.5625
        0.78452 0.61887 0.875
        0.68097 0.0 0.5
        0.18097 0.5 0.5
        0.22957 0.08655 0.9375
        0.67463 0.36887 0.75
        0.39322 0.37776 0.75
        0.05195 0.46844 0.75
        0.44805 0.19164 0.9375
        0.44465 0.43889 0.0625
        0.27043 0.58655 0.9375
        0.35416 0.25966 0.0625
        0.64243 0.10234 0.125
        """

        self.coord = "relative"

        self.cages = """
        14 -0.1333 -0.37072 0.75
        14 -0.1333 0.37072 0.25
        12 0.0 0.5 0.0
        14 0.6333 0.87072 0.25
        12 0.14692 0.0 0.0
        15 0.77136 0.34618 -0.25
        14 0.3667 0.87072 0.25
        14 0.1333 -0.37072 -0.25
        15 -0.27136 -0.15382 0.75
        16 0.5 0.61825 0.25
        12 0.5 0.0 0.5
        15 -0.27136 0.15382 0.25
        12 0.64692 0.5 -0.5
        12 0.42551 -0.25554 -0.25
        15 0.22864 0.34618 0.75
        12 -0.14692 0.0 0.0
        12 0.92551 0.75554 0.25
        12 0.07449 0.75554 0.25
        12 0.07449 0.24446 0.75
        12 0.35308 0.5 0.5
        14 0.3667 0.12928 0.75
        12 0.5 1.0 0.0
        14 0.1333 0.37072 0.25
        12 0.64692 0.5 0.0
        15 0.22864 0.65382 0.25
        12 0.92551 0.24446 -0.25
        12 0.35308 0.5 0.0
        14 0.6333 0.12928 -0.25
        15 0.27136 -0.15382 -0.25
        16 0.0 0.11825 0.25
        12 0.0 -0.5 0.5
        16 0.0 -0.11825 0.75
        15 0.77136 0.65382 0.25
        12 0.14692 0.0 -0.5
        16 0.5 0.38175 0.75
        12 -0.42551 0.25554 0.25
        12 -0.14692 0.0 0.5
        12 0.42551 0.25554 0.25
        12 -0.42551 -0.25554 0.75
        15 0.27136 0.15382 0.25
        """

        self.bondlen = 3

        self.cell = """
        49.88351534567095 17.84622037904867 13.527667556276105
        """

        self.density = 0.565900889084785

        self.cell = cellvectors(a=49.88351534567095,
                                b=17.84622037904867,
                                c=13.527667556276105)
