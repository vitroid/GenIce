# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (20,10,4,6,)
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
        143 214
        126 187
        15 56
        103 200
        185 116
        140 93
        24 195
        38 18
        36 30
        88 11
        226 212
        110 63
        102 209
        115 200
        155 216
        33 40
        62 90
        128 207
        122 206
        186 57
        49 74
        89 69
        62 58
        104 207
        36 207
        12 75
        23 189
        211 180
        3 136
        116 10
        0 197
        139 223
        65 53
        129 194
        67 153
        135 77
        26 2
        101 124
        133 144
        54 108
        13 158
        183 41
        12 159
        20 10
        224 56
        146 49
        78 200
        172 136
        74 145
        221 199
        9 22
        36 50
        49 46
        177 227
        79 213
        64 50
        5 101
        110 100
        113 193
        211 70
        86 57
        118 59
        134 78
        14 178
        50 208
        161 191
        115 149
        173 28
        15 62
        140 153
        143 45
        27 182
        167 201
        203 130
        17 149
        128 161
        37 179
        88 66
        189 176
        132 4
        45 115
        118 125
        21 167
        166 164
        156 226
        1 0
        17 16
        68 215
        87 86
        62 19
        204 84
        142 109
        45 135
        66 219
        35 71
        13 100
        24 141
        63 97
        180 191
        79 91
        119 55
        154 114
        146 172
        80 125
        122 33
        89 79
        172 185
        18 4
        196 181
        126 93
        105 148
        193 199
        30 80
        169 113
        9 108
        70 163
        6 197
        44 202
        138 42
        160 178
        74 193
        204 132
        169 136
        174 215
        54 21
        9 94
        181 97
        197 175
        69 66
        138 29
        65 78
        135 127
        184 225
        49 96
        35 162
        115 120
        25 98
        137 6
        28 214
        16 119
        73 224
        210 165
        137 147
        224 217
        13 202
        222 68
        78 215
        85 8
        3 50
        133 85
        161 76
        138 175
        139 175
        188 124
        178 124
        142 113
        1 148
        44 176
        22 40
        35 183
        37 107
        80 51
        134 28
        174 206
        189 141
        152 166
        103 174
        36 191
        107 51
        212 225
        226 75
        221 173
        156 81
        75 176
        1 13
        114 117
        55 186
        222 19
        61 68
        12 23
        25 209
        87 114
        156 159
        223 100
        95 56
        104 164
        77 165
        87 34
        184 166
        187 213
        105 74
        189 217
        30 151
        195 97
        132 107
        60 218
        127 222
        3 216
        169 145
        129 46
        70 83
        16 34
        48 24
        44 171
        26 17
        52 109
        118 166
        29 131
        126 120
        92 87
        121 131
        198 73
        20 16
        77 227
        142 219
        23 58
        157 41
        117 181
        205 159
        113 185
        60 187
        54 162
        2 46
        160 118
        193 123
        182 84
        203 83
        137 91
        154 57
        19 122
        55 199
        168 192
        106 33
        220 98
        96 191
        95 192
        190 101
        31 42
        151 145
        81 164
        129 150
        77 218
        125 207
        129 163
        134 120
        99 190
        99 192
        76 225
        54 40
        94 206
        24 121
        108 71
        138 1
        107 64
        119 86
        170 184
        144 212
        39 186
        71 112
        7 22
        94 167
        34 102
        14 47
        64 152
        70 140
        157 187
        73 121
        119 173
        136 31
        90 220
        37 208
        38 209
        155 211
        25 188
        11 116
        133 75
        34 98
        149 214
        179 42
        188 90
        104 64
        69 153
        211 76
        27 32
        44 110
        128 59
        17 203
        14 32
        91 158
        38 204
        209 57
        45 111
        96 151
        106 18
        25 53
        84 51
        160 8
        168 131
        127 41
        15 47
        225 81
        175 212
        218 111
        61 206
        29 226
        92 32
        6 76
        213 150
        43 100
        190 32
        3 180
        168 85
        154 196
        117 4
        35 109
        221 72
        95 14
        7 196
        208 82
        10 112
        126 111
        220 200
        146 83
        108 186
        124 204
        202 131
        73 192
        63 51
        162 165
        178 182
        39 102
        5 132
        183 227
        210 10
        67 213
        190 90
        170 147
        94 177
        99 58
        0 88
        26 120
        171 217
        221 46
        198 122
        7 106
        23 56
        205 101
        155 194
        171 121
        103 19
        92 220
        67 163
        198 58
        43 148
        0 91
        47 65
        72 214
        135 103
        20 130
        174 227
        59 170
        142 89
        61 224
        43 145
        2 83
        205 164
        188 47
        99 159
        179 110
        82 43
        105 11
        27 117
        218 66
        88 52
        21 210
        33 167
        201 114
        93 150
        2 93
        37 195
        38 196
        60 52
        125 182
        216 151
        168 156
        128 81
        139 158
        98 149
        105 172
        20 39
        21 39
        95 8
        155 184
        137 150
        79 52
        60 183
        203 72
        48 179
        130 199
        29 48
        133 202
        12 8
        161 147
        216 152
        96 194
        147 194
        127 162
        208 31
        170 144
        6 140
        55 112
        7 217
        219 165
        157 134
        222 40
        146 180
        41 215
        15 68
        106 141
        59 85
        123 112
        72 163
        22 61
        63 82
        80 152
        160 205
        109 123
        201 102
        48 176
        198 141
        30 82
        26 173
        86 53
        223 42
        111 153
        4 195
        157 143
        71 177
        69 197
        11 123
        18 201
        31 148
        158 144
        9 154
        210 177
        104 5
        53 28
        171 97
        92 65
        27 5
        185 130
        143 67
        139 89
        84 181
        223 169
        219 116
        """

        self.waters = """
        0.5 0.19693 0.86923
        0.5 0.26643 0.94069
        0.5 0.77187 0.73445
        0.81248 0.61786 0.93874
        0.68454 0.51306 0.25765
        0.68454 0.70923 0.2135
        0.625 0.98985 0.87229
        0.125 0.35558 0.32428
        0.19395 0.97882 0.18494
        0.31248 0.38214 0.43874
        0.625 0.44423 0.63688
        0.5 0.3763 0.77139
        0.0 0.04911 0.18723
        0.30954 0.24325 0.98383
        0.31546 0.86945 0.29698
        0.125 0.01015 0.37229
        0.625 0.62074 0.58038
        0.69395 0.70958 0.62413
        0.80954 0.4799 0.32581
        0.81253 0.08507 0.42581
        0.75 0.53239 0.61175
        0.81253 0.39622 0.52165
        0.18748 0.30419 0.39623
        0.0 0.10182 0.25288
        0.68105 0.32408 0.19759
        0.0 0.73357 0.44069
        0.5 0.76047 0.65413
        0.5 0.7024 0.25864
        0.18105 0.79924 0.56856
        0.68105 0.20077 0.06856
        0.31253 0.60378 0.02165
        0.69046 0.43302 0.96507
        0.5 0.80836 0.31371
        0.81253 0.30419 0.39623
        0.68454 0.64788 0.49053
        0.18748 0.26651 0.62529
        0.5 0.67141 0.01747
        0.75 0.46761 0.11175
        0.0 0.55915 0.33421
        0.875 0.50181 0.53958
        0.0 0.26945 0.44142
        0.18748 0.07805 0.56302
        0.81546 0.35212 0.99053
        0.30954 0.43302 0.96507
        0.19395 0.29042 0.12413
        0.81895 0.98452 0.61095
        0.30605 0.71891 0.76633
        0.19046 0.88447 0.37382
        0.80605 0.29042 0.12413
        0.375 0.64442 0.82428
        0.68748 0.60378 0.02165
        0.125 0.55577 0.13688
        0.30954 0.20512 0.74725
        0.19046 0.75676 0.48383
        0.0 0.32859 0.51747
        0.25 0.53239 0.61175
        0.19395 0.06542 0.29774
        0.19046 0.56698 0.46507
        0.80605 0.06542 0.29774
        0.31253 0.92195 0.06302
        0.375 0.13214 0.68008
        0.31248 0.1933 0.38314
        0.875 0.01015 0.37229
        0.25 0.46761 0.11175
        0.81248 0.64262 0.08383
        0.31546 0.84656 0.43369
        0.69046 0.20512 0.74725
        0.0 0.95607 0.73317
        0.18748 0.08507 0.42581
        0.81546 0.13056 0.79698
        0.81895 0.82492 0.80607
        0.31248 0.35738 0.58383
        0.0 0.74554 0.69247
        0.5 0.16719 0.26733
        0.30954 0.5201 0.82581
        0.0 0.12511 0.12211
        0.68748 0.91493 0.92581
        0.68753 0.16897 0.6086
        0.375 0.93385 0.49152
        0.18454 0.13056 0.79698
        0.18753 0.64262 0.08383
        0.68748 0.92195 0.06302
        0.375 0.49819 0.03958
        0.69395 0.71891 0.76633
        0.19046 0.59914 0.22091
        0.31895 0.01548 0.11095
        0.31546 0.64788 0.49053
        0.5 0.65629 0.4358
        0.5 0.25262 0.78981
        0.0 0.19165 0.81371
        0.80954 0.88447 0.37382
        0.30954 0.11554 0.87382
        0.5 0.78404 0.39942
        0.5 0.89819 0.75288
        0.5 0.32587 0.45785
        0.31895 0.98724 0.25671
        0.31253 0.69582 0.89623
        0.375 0.43173 0.18738
        0.80954 0.75676 0.48383
        0.68105 0.98724 0.25671
        0.18454 0.35212 0.99053
        0.80954 0.79488 0.24725
        0.80954 0.56698 0.46507
        0.68753 0.04668 0.48517
        0.68748 0.73349 0.12529
        0.5 0.44085 0.83421
        0.875 0.35558 0.32428
        0.875 0.55577 0.13688
        0.18748 0.39622 0.52165
        0.18454 0.29077 0.7135
        0.125 0.37926 0.08038
        0.69395 0.02119 0.68494
        0.375 0.44423 0.63688
        0.0 0.42435 0.7928
        0.5 0.55385 0.38426
        0.69395 0.89493 0.57166
        0.69046 0.40086 0.72091
        0.5 0.57565 0.2928
        0.18753 0.83103 0.1086
        0.375 0.62074 0.58038
        0.5 0.8749 0.62211
        0.5 0.25446 0.19247
        0.68753 0.1933 0.38314
        0.30954 0.40086 0.72091
        0.0 0.74738 0.28981
        0.31253 0.73349 0.12529
        0.5 0.95089 0.68723
        0.0 0.12677 0.54199
        0.5 0.87323 0.04199
        0.18105 0.82492 0.80607
        0.875 0.56827 0.68738
        0.5 0.18147 0.11014
        0.80954 0.59914 0.22091
        0.19395 0.10507 0.07166
        0.30605 0.89493 0.57166
        0.81253 0.07805 0.56302
        0.81546 0.50387 0.90194
        0.375 0.98985 0.87229
        0.69046 0.24325 0.98383
        0.0 0.21596 0.89942
        0.69395 0.93458 0.79774
        0.80605 0.28109 0.26633
        0.0 0.2976 0.75864
        0.0 0.93971 0.63785
        0.125 0.06615 0.99152
        0.18454 0.50387 0.90194
        0.625 0.64442 0.82428
        0.31253 0.91493 0.92581
        0.5 0.38666 0.92409
        0.81895 0.79924 0.56856
        0.30605 0.93458 0.79774
        0.18753 0.61786 0.93874
        0.0 0.69144 0.04611
        0.81895 0.01277 0.75671
        0.31546 0.49613 0.40194
        0.0 0.79638 0.92146
        0.68105 0.01548 0.11095
        0.18105 0.98452 0.61095
        0.18454 0.15344 0.93369
        0.80605 0.97882 0.18494
        0.125 0.86786 0.18008
        0.5 0.86032 0.95802
        0.0 0.24908 0.58223
        0.0 0.83282 0.76733
        0.81248 0.83103 0.1086
        0.81253 0.26651 0.62529
        0.0 0.80779 0.06296
        0.68753 0.38214 0.43874
        0.5 0.06029 0.13785
        0.0 0.44615 0.88426
        0.18753 0.95332 0.98517
        0.31895 0.32408 0.19759
        0.69046 0.5201 0.82581
        0.30605 0.70958 0.62413
        0.5 0.12017 0.48769
        0.81546 0.15344 0.93369
        0.0 0.23954 0.15413
        0.5 0.30856 0.54611
        0.19046 0.79488 0.24725
        0.875 0.37926 0.08038
        0.68748 0.69582 0.89623
        0.31546 0.51306 0.25765
        0.31546 0.70923 0.2135
        0.31248 0.16897 0.6086
        0.0 0.87983 0.98769
        0.81546 0.48694 0.75765
        0.125 0.50181 0.53958
        0.30605 0.02119 0.68494
        0.0 0.80307 0.36923
        0.0 0.22813 0.23445
        0.68454 0.86945 0.29698
        0.5 0.73055 0.94142
        0.5 0.04394 0.23317
        0.18454 0.48694 0.75765
        0.18753 0.8067 0.88314
        0.625 0.43173 0.18738
        0.19046 0.4799 0.32581
        0.69046 0.11554 0.87382
        0.68105 0.17508 0.30607
        0.125 0.56827 0.68738
        0.625 0.93385 0.49152
        0.68454 0.49613 0.40194
        0.31895 0.20077 0.06856
        0.81895 0.67592 0.69759
        0.0 0.6237 0.27139
        0.875 0.86786 0.18008
        0.5 0.20363 0.42146
        0.5 0.75092 0.08223
        0.625 0.49819 0.03958
        0.0 0.61334 0.42409
        0.68753 0.35738 0.58383
        0.81248 0.8067 0.88314
        0.875 0.06615 0.99152
        0.18105 0.01277 0.75671
        0.0 0.81853 0.61014
        0.31248 0.04668 0.48517
        0.0 0.67413 0.95785
        0.19395 0.28109 0.26633
        0.625 0.13214 0.68008
        0.81546 0.29077 0.7135
        0.68454 0.84656 0.43369
        0.18105 0.67592 0.69759
        0.0 0.13968 0.45802
        0.0 0.34372 0.9358
        0.31895 0.17508 0.30607
        0.81248 0.95332 0.98517
        0.80605 0.10507 0.07166
        0.5 0.19221 0.56296
        """

        self.coord = "relative"

        self.cages = """
        12 0.0 -0.02863 1.08298
        14 0.0 0.41877 0.64974
        16 0.5 -0.92953 0.98878
        12 0.0 0.22781 0.02427
        12 0.2242 -0.15062 0.69059
        12 0.2499 0.2235 0.49886
        12 0.0 -0.07618 1.28476
        12 0.7499 -0.2235 0.99886
        16 0.0 -0.00505 0.87004
        12 0.7242 0.15062 1.19059
        14 0.5 0.35416 0.0706
        15 0.5 -0.42909 0.71107
        15 0.0 0.42909 0.21107
        14 0.26184 0.67946 0.35168
        12 0.5 0.02863 0.58298
        14 0.0 0.50597 0.01558
        12 0.5 -0.44859 0.92241
        12 0.5 -0.1822 0.84556
        14 0.76184 -0.67946 0.85168
        16 0.5 0.00505 0.37004
        12 0.2501 -0.2235 0.99886
        12 0.0 0.44859 0.42241
        12 0.2758 0.15062 0.19059
        12 0.7501 0.2235 1.49886
        15 0.5 -0.13055 1.18144
        16 0.0 -0.36245 0.81809
        15 0.0 0.13055 0.68144
        14 0.73816 0.67946 1.35168
        16 0.5 0.36245 0.31809
        12 0.5 -0.70679 0.67113
        12 0.0 0.70679 0.17113
        16 0.0 0.92953 0.48878
        12 0.5 0.07618 0.78476
        14 0.5 -0.41877 1.14974
        12 0.5 -0.22781 0.52427
        12 0.0 0.1822 0.34556
        12 0.7758 -0.15062 0.69059
        14 0.23816 -0.67946 0.85168
        14 0.5 -0.50597 0.51558
        14 0.0 -0.35416 0.5706
        """

        self.bondlen = 3

        self.cell = """
        14.082693759637738 25.66558501658822 36.64668377743608
        """

        self.density = 0.5145114179843194

        self.cell = cellvectors(a=14.082693759637738,
                                b=25.66558501658822,
                                c=36.64668377743608)