# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (22,2,8,6,)
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
        104 202
        1 71
        214 85
        129 111
        103 55
        16 71
        108 181
        30 22
        86 114
        133 137
        69 191
        204 61
        1 95
        26 4
        124 147
        123 198
        151 70
        99 150
        25 114
        55 48
        188 197
        16 96
        162 129
        112 66
        133 163
        91 119
        119 92
        85 164
        157 110
        172 197
        124 19
        154 181
        121 83
        150 153
        35 137
        1 169
        187 162
        82 47
        24 156
        28 178
        108 102
        160 137
        144 215
        93 110
        162 177
        122 139
        82 59
        131 198
        74 85
        182 129
        25 215
        10 13
        122 42
        187 69
        213 100
        157 199
        42 155
        174 77
        5 45
        11 54
        21 134
        5 3
        175 19
        91 28
        77 72
        152 177
        26 86
        186 176
        45 34
        30 189
        139 118
        144 178
        130 113
        21 173
        135 170
        161 61
        185 146
        21 66
        159 87
        121 78
        192 146
        167 170
        122 138
        98 84
        19 153
        209 152
        0 102
        2 103
        188 173
        97 115
        213 196
        4 84
        47 39
        151 175
        169 109
        44 132
        201 23
        90 149
        27 168
        48 3
        7 28
        187 41
        92 158
        122 188
        51 160
        105 156
        100 212
        148 18
        72 38
        112 152
        139 174
        154 89
        181 164
        128 29
        28 25
        43 67
        188 94
        50 68
        12 32
        31 193
        181 94
        210 6
        81 156
        0 93
        97 153
        60 14
        201 198
        186 53
        185 53
        144 41
        189 71
        189 72
        195 13
        99 103
        57 196
        21 138
        175 103
        11 163
        199 59
        194 69
        156 202
        57 162
        176 54
        101 25
        136 212
        113 115
        144 194
        109 37
        127 202
        142 150
        59 178
        141 56
        118 183
        6 106
        206 64
        125 149
        27 74
        115 100
        199 104
        75 68
        195 111
        112 65
        36 58
        107 179
        51 15
        109 23
        125 141
        0 128
        46 179
        192 207
        64 32
        109 116
        40 182
        182 41
        16 140
        45 113
        63 167
        40 215
        88 54
        57 8
        56 8
        128 111
        89 43
        80 52
        22 120
        193 208
        140 95
        206 201
        36 193
        65 106
        64 196
        27 3
        18 60
        203 171
        191 68
        97 147
        12 191
        130 46
        130 44
        183 140
        2 126
        161 86
        73 75
        30 126
        159 143
        120 14
        163 165
        1 33
        36 171
        121 211
        191 131
        134 24
        140 9
        7 82
        157 91
        93 105
        78 167
        143 11
        155 114
        20 213
        74 55
        136 98
        38 92
        149 129
        200 166
        17 119
        185 17
        185 18
        51 94
        67 178
        138 61
        22 9
        49 48
        54 39
        210 203
        7 155
        117 172
        29 173
        87 146
        43 91
        0 182
        37 6
        146 76
        127 160
        211 172
        84 152
        133 76
        176 38
        20 107
        136 184
        186 180
        80 105
        83 47
        148 120
        118 85
        209 10
        11 104
        176 17
        78 197
        107 64
        49 99
        34 70
        153 212
        159 165
        89 102
        49 97
        5 46
        73 63
        174 211
        65 170
        184 150
        131 40
        46 76
        123 86
        159 154
        127 163
        154 160
        31 190
        141 12
        184 10
        88 180
        70 99
        40 79
        35 147
        74 124
        167 114
        171 126
        78 42
        87 168
        190 53
        62 55
        95 19
        90 195
        214 62
        61 155
        24 108
        87 137
        161 66
        70 10
        121 204
        190 107
        88 18
        81 67
        151 116
        110 202
        30 58
        143 17
        82 204
        89 110
        77 83
        200 23
        164 165
        45 49
        179 180
        57 12
        27 165
        184 210
        36 132
        143 43
        84 111
        125 206
        209 116
        168 147
        128 80
        31 44
        210 200
        207 115
        205 194
        157 101
        4 29
        209 6
        130 192
        148 58
        83 14
        158 14
        101 79
        133 3
        62 15
        124 183
        38 47
        73 65
        193 33
        2 9
        141 44
        31 206
        16 120
        108 173
        139 96
        50 32
        32 132
        135 52
        63 205
        205 59
        51 117
        118 94
        63 204
        134 52
        8 13
        145 127
        50 37
        80 177
        168 207
        161 73
        187 105
        35 15
        20 207
        169 126
        13 177
        15 183
        136 195
        58 180
        135 197
        39 60
        90 100
        96 72
        56 113
        166 106
        50 201
        125 131
        149 196
        119 60
        5 20
        145 117
        174 214
        88 76
        7 92
        35 48
        211 138
        26 170
        142 203
        200 98
        169 175
        203 33
        123 215
        52 112
        62 9
        22 77
        135 29
        4 66
        148 208
        123 75
        208 53
        199 39
        96 158
        69 79
        34 8
        71 208
        214 117
        198 166
        151 212
        93 79
        102 81
        158 42
        145 24
        205 101
        142 95
        90 56
        142 2
        23 33
        81 41
        34 213
        145 164
        68 106
        104 67
        190 192
        186 189
        132 179
        166 26
        171 37
        134 172
        98 116
        75 194
        """

        self.waters="""
        0.85656 0.93312 0.86012
        0.40941 0.06689 0.29978
        0.5011 0.55812 0.32846
        0.84593 0.43763 0.40153
        0.55993 0.05812 0.92548
        0.90586 0.43312 0.31968
        0.44007 0.55812 0.07452
        0.25 0.125 0.69948
        0.85141 0.55812 0.105
        0.49866 0.625 0.42519
        0.64859 0.55812 0.105
        0.99866 0.375 0.57482
        0.05993 0.55812 0.07452
        0.75 0.625 0.05449
        0.31263 0.75 0.56016
        0.65432 0.75 0.49443
        0.40696 0.93763 0.45457
        0.0917 0.06263 0.52938
        0.15541 0.75 0.46944
        0.59414 0.06689 0.31968
        0.93248 0.25 0.27614
        0.58789 0.25 0.79383
        0.40696 0.56237 0.45457
        0.35656 0.06689 0.13988
        0.75 0.375 0.74256
        0.19007 0.93312 0.81686
        0.43493 0.93312 0.93508
        0.81263 0.25 0.43985
        0.15677 0.05812 0.74302
        0.64345 0.93312 0.86012
        0.34434 0.43737 0.37269
        0.15297 0.94189 0.22102
        0.14345 0.43312 0.13988
        0.34704 0.94189 0.22102
        0.80993 0.43312 0.18315
        0.75 0.75 0.44124
        0.25 0.625 0.25744
        0.35656 0.43312 0.13988
        0.25 0.25 0.55877
        0.15407 0.56237 0.59847
        0.06507 0.93312 0.93508
        0.96366 0.25 0.88077
        0.40586 0.93312 0.68033
        0.0011 0.05812 0.67154
        0.08789 0.75 0.20617
        0.84323 0.55812 0.25699
        0.9989 0.55812 0.32846
        0.25 0.43737 0.61298
        0.75 0.56263 0.38702
        0.75 0.625 0.30052
        0.25 0.375 0.09516
        0.68763 0.75 0.59081
        0.64345 0.56689 0.86012
        0.15567 0.06263 0.37269
        0.0917 0.43737 0.52938
        0.65407 0.43763 0.40153
        0.91134 0.75 0.1353
        0.93493 0.43312 0.06493
        0.25 0.56237 0.35107
        0.15677 0.44189 0.74302
        0.18737 0.75 0.56016
        0.40316 0.25 0.77065
        0.5917 0.56263 0.47062
        0.30993 0.56689 0.81686
        0.07837 0.25 0.16059
        0.43493 0.56689 0.93508
        0.53634 0.25 0.88077
        0.0038 0.25 0.72959
        0.25 0.5 0.01539
        0.06507 0.56689 0.93508
        0.69007 0.43312 0.18315
        0.34434 0.06263 0.37269
        0.34568 0.25 0.50558
        0.35141 0.44189 0.895
        0.68737 0.25 0.43985
        0.25 0.375 0.94552
        0.00134 0.625 0.42519
        0.40831 0.43737 0.52938
        0.43248 0.75 0.72386
        0.04663 0.75 0.88456
        0.75 0.625 0.90484
        0.91211 0.25 0.79383
        0.25 0.375 0.69948
        0.34593 0.56237 0.59847
        0.625 0.125 0.99041
        0.65541 0.25 0.53056
        0.35141 0.05812 0.895
        0.90831 0.93737 0.47062
        0.09304 0.56237 0.45457
        0.90941 0.93312 0.70023
        0.85141 0.94189 0.105
        0.09414 0.93312 0.68033
        0.25 0.06263 0.61298
        0.92163 0.75 0.83941
        0.65567 0.93737 0.62731
        0.5011 0.94189 0.32846
        0.40831 0.06263 0.52938
        0.75 0.875 0.30052
        0.56507 0.06689 0.06493
        0.65677 0.55812 0.25699
        0.80993 0.06689 0.18315
        0.13015 0.75 0.81221
        0.84704 0.05812 0.77898
        0.59414 0.43312 0.31968
        0.0011 0.44189 0.67154
        0.85656 0.56689 0.86012
        0.375 0.625 0.00959
        0.05748 0.25 0.26336
        0.75 0.125 0.74256
        0.42163 0.25 0.16059
        0.94252 0.75 0.73665
        0.75 0.0 0.98461
        0.55993 0.44189 0.92548
        0.90316 0.75 0.22935
        0.30993 0.93312 0.81686
        0.84323 0.94189 0.25699
        0.54663 0.25 0.11544
        0.65567 0.56263 0.62731
        0.59304 0.06237 0.54543
        0.15407 0.93763 0.59847
        0.34459 0.75 0.46944
        0.40586 0.56689 0.68033
        0.4989 0.05812 0.67154
        0.25 0.125 0.94552
        0.65407 0.06237 0.40153
        0.05993 0.94189 0.07452
        0.40941 0.43312 0.29978
        0.84434 0.56263 0.62731
        0.75 0.875 0.90484
        0.875 0.125 0.99041
        0.9962 0.75 0.27042
        0.125 0.875 0.00959
        0.15297 0.55812 0.22102
        0.90831 0.56263 0.47062
        0.65297 0.44189 0.77898
        0.57837 0.75 0.83941
        0.64859 0.94189 0.105
        0.84568 0.75 0.49443
        0.4962 0.25 0.72959
        0.50134 0.125 0.57482
        0.49866 0.875 0.42519
        0.03634 0.75 0.11924
        0.5038 0.75 0.27042
        0.99866 0.125 0.57482
        0.08866 0.25 0.8647
        0.75 0.43763 0.64893
        0.00134 0.875 0.42519
        0.75 0.93737 0.38702
        0.25 0.75 0.40914
        0.93493 0.06689 0.06493
        0.59684 0.75 0.22935
        0.63015 0.25 0.18779
        0.625 0.375 0.99041
        0.65677 0.94189 0.25699
        0.84434 0.93737 0.62731
        0.34323 0.05812 0.74302
        0.84704 0.44189 0.77898
        0.06752 0.75 0.72386
        0.34593 0.93763 0.59847
        0.90696 0.06237 0.54543
        0.81237 0.75 0.59081
        0.41134 0.25 0.8647
        0.875 0.375 0.99041
        0.90696 0.43763 0.54543
        0.75 0.25 0.59087
        0.84459 0.25 0.53056
        0.375 0.875 0.00959
        0.36986 0.75 0.81221
        0.84593 0.06237 0.40153
        0.44252 0.25 0.26336
        0.45337 0.75 0.88456
        0.34704 0.55812 0.22102
        0.59059 0.56689 0.70023
        0.65297 0.05812 0.77898
        0.50134 0.375 0.57482
        0.56752 0.25 0.27614
        0.15432 0.25 0.50558
        0.75 0.5 0.98461
        0.09684 0.25 0.77065
        0.09059 0.43312 0.29978
        0.15567 0.43737 0.37269
        0.75 0.06237 0.64893
        0.94007 0.05812 0.92548
        0.5917 0.93737 0.47062
        0.58866 0.75 0.1353
        0.09304 0.93763 0.45457
        0.18763 0.25 0.40919
        0.94007 0.44189 0.92548
        0.59059 0.93312 0.70023
        0.31237 0.25 0.40919
        0.09059 0.06689 0.29978
        0.125 0.625 0.00959
        0.9989 0.94189 0.32846
        0.25 0.875 0.25744
        0.14859 0.44189 0.895
        0.75 0.875 0.05449
        0.95337 0.25 0.11544
        0.55748 0.75 0.73665
        0.25 0.0 0.01539
        0.09414 0.56689 0.68033
        0.44007 0.94189 0.07452
        0.25 0.125 0.09516
        0.90941 0.56689 0.70023
        0.41211 0.75 0.20617
        0.34323 0.44189 0.74302
        0.19007 0.56689 0.81686
        0.14345 0.06689 0.13988
        0.90586 0.06689 0.31968
        0.25 0.93763 0.35107
        0.56507 0.43312 0.06493
        0.46366 0.75 0.11924
        0.4989 0.44189 0.67154
        0.69007 0.06689 0.18315
        0.86986 0.25 0.18779
        0.59304 0.43763 0.54543
        0.14859 0.05812 0.895
        """

        self.coord= "relative"

        self.cages="""
        16 -0.00051 0.75 -0.42264
        12 0.25 0.00053 0.48213
        15 -0.08407 0.75 -0.01159
        12 0.01029 0.98246 -0.19046
        12 0.38322 0.25 0.6354
        12 -0.01029 0.48246 0.19046
        12 -0.25 0.50053 -0.48213
        12 -0.37785 0.25 -0.33614
        15 -0.41593 0.75 -0.01159
        12 0.37785 -0.25 0.33614
        15 -0.75 0.75 -0.301
        12 -0.11678 0.75 -0.6354
        16 0.00051 0.25 0.42264
        12 -0.38322 0.75 -0.6354
        15 0.08407 0.25 0.01159
        12 0.12215 0.75 0.33614
        16 -0.25 0.75 -0.24986
        12 -0.51029 0.51754 -0.19046
        12 0.25 0.49947 0.48213
        12 0.51029 -0.51754 0.19046
        16 0.25 0.25 0.24986
        14 -0.75 0.75 -0.06924
        16 -0.49949 0.75 -0.42264
        12 0.11678 0.25 0.6354
        15 0.25 0.75 0.10761
        15 0.75 0.25 0.301
        12 -0.25 -0.00053 -0.48213
        14 0.75 0.25 0.06924
        12 0.01029 -0.48246 -0.19046
        12 -0.01029 -0.98246 0.19046
        15 0.41593 0.25 0.01159
        15 -0.25 1.25 -0.10761
        16 0.49949 0.25 0.42264
        12 0.25 0.25 0.82812
        12 0.51029 0.01754 0.19046
        12 -0.12215 1.25 -0.33614
        12 -0.25 0.75 -0.82812
        12 -0.51029 -0.01754 -0.19046
        """

        self.bondlen = 3


        self.cell = """
        26.304377127143944 14.227456687621972 34.12719015977335
        """

        self.density = 0.5055099148935727



        self.cell = cellvectors(a=26.304377127143944,
                           b=14.227456687621972,
                           c=34.12719015977335)
