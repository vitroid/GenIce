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
        157 40
        141 41
        188 88
        75 85
        115 210
        183 82
        177 87
        62 211
        95 48
        189 43
        143 80
        110 76
        226 213
        56 48
        70 85
        193 162
        114 20
        30 194
        24 16
        7 79
        113 96
        166 136
        226 187
        78 175
        127 60
        142 215
        3 59
        3 58
        136 22
        119 52
        147 44
        188 28
        26 54
        183 20
        112 189
        226 181
        64 131
        45 201
        10 154
        22 116
        149 70
        2 42
        27 196
        6 11
        118 128
        165 30
        174 216
        118 142
        183 98
        66 103
        13 167
        126 5
        113 11
        0 11
        7 170
        161 162
        91 16
        216 138
        123 25
        76 225
        41 192
        206 122
        123 82
        177 34
        58 207
        140 102
        94 30
        124 69
        124 66
        186 152
        196 210
        174 69
        44 176
        81 94
        127 86
        32 134
        222 41
        166 36
        38 141
        29 87
        3 107
        17 80
        54 190
        62 180
        65 182
        147 83
        99 100
        10 49
        124 138
        177 179
        8 113
        185 216
        160 225
        13 155
        154 97
        199 5
        203 105
        204 197
        204 199
        32 217
        204 92
        59 1
        125 178
        150 40
        205 96
        75 166
        19 152
        58 128
        77 208
        26 51
        38 104
        172 131
        220 5
        135 164
        209 57
        150 9
        90 61
        47 136
        196 53
        177 133
        126 48
        138 176
        209 64
        227 158
        223 221
        198 74
        116 12
        115 12
        117 167
        171 111
        23 171
        73 91
        132 163
        26 169
        85 137
        106 13
        193 65
        210 22
        195 180
        72 132
        72 131
        152 194
        69 211
        122 94
        70 145
        77 4
        169 39
        140 71
        121 132
        224 95
        13 181
        77 191
        164 162
        97 56
        15 187
        218 175
        46 22
        192 51
        107 178
        67 111
        156 86
        143 103
        144 91
        2 188
        148 124
        1 159
        155 185
        106 18
        126 151
        189 69
        170 212
        38 190
        112 50
        31 163
        18 117
        50 156
        102 180
        112 127
        106 110
        31 21
        28 82
        71 200
        128 9
        36 217
        6 47
        171 115
        36 191
        27 47
        15 207
        112 17
        107 201
        174 74
        214 109
        215 110
        224 94
        108 147
        185 109
        184 144
        50 92
        132 79
        123 184
        42 78
        160 32
        111 104
        18 174
        15 157
        68 56
        43 60
        87 199
        12 104
        46 227
        192 146
        137 28
        108 90
        134 191
        125 37
        163 221
        42 86
        155 203
        159 208
        148 156
        115 4
        40 163
        27 150
        88 212
        99 68
        24 206
        29 131
        101 61
        183 196
        10 144
        186 92
        215 52
        119 203
        212 111
        206 101
        129 191
        172 55
        217 104
        46 63
        133 92
        19 179
        63 128
        19 204
        14 213
        160 102
        134 41
        179 81
        84 181
        46 53
        173 168
        219 193
        118 120
        67 221
        176 17
        57 212
        34 197
        37 201
        205 224
        192 190
        149 60
        155 169
        113 82
        87 168
        58 52
        8 55
        32 146
        6 55
        161 89
        151 205
        222 198
        142 225
        33 122
        39 198
        34 202
        125 187
        0 25
        62 89
        98 21
        0 99
        142 146
        135 122
        103 211
        35 52
        149 93
        148 91
        101 97
        116 1
        114 172
        107 129
        141 208
        209 139
        90 213
        222 71
        89 219
        171 188
        83 101
        33 61
        89 222
        36 116
        19 126
        117 81
        63 157
        64 158
        65 159
        202 175
        99 96
        219 102
        61 194
        2 218
        130 199
        221 57
        66 44
        8 98
        21 158
        6 20
        223 157
        125 159
        220 98
        137 11
        71 195
        170 217
        214 81
        220 153
        145 73
        1 129
        224 68
        49 127
        137 136
        186 185
        140 141
        121 120
        49 80
        201 162
        120 57
        130 153
        62 74
        151 97
        10 83
        24 100
        68 206
        27 4
        18 109
        195 176
        106 35
        154 96
        145 184
        108 164
        60 73
        114 220
        2 184
        186 167
        25 20
        139 130
        93 88
        100 73
        72 175
        14 105
        143 164
        173 205
        189 133
        135 165
        7 146
        25 56
        51 198
        7 120
        40 53
        55 158
        30 167
        35 105
        200 182
        31 153
        100 49
        150 207
        33 214
        67 190
        156 197
        59 187
        144 86
        95 168
        54 225
        66 16
        83 16
        147 80
        151 194
        85 88
        23 218
        117 133
        200 44
        84 35
        216 39
        64 121
        31 209
        213 37
        9 227
        0 145
        148 43
        208 12
        26 119
        154 123
        54 223
        15 119
        226 203
        110 74
        70 218
        214 152
        14 45
        135 45
        202 29
        93 139
        165 181
        50 138
        8 173
        169 76
        3 84
        53 21
        207 178
        149 202
        134 219
        193 129
        37 182
        42 93
        23 75
        78 130
        47 227
        75 170
        108 182
        211 17
        105 109
        72 153
        29 139
        118 223
        140 65
        77 178
        23 79
        43 34
        172 168
        114 48
        210 28
        59 63
        195 39
        38 160
        197 78
        9 121
        173 5
        4 166
        165 90
        161 103
        215 51
        67 79
        84 45
        95 179
        14 33
        180 76
        143 24
        161 200
        """

        self.waters = """
        0.18869 0.72845 0.81382
        0.5 0.33257 0.88262
        0.81685 0.67268 0.61855
        0.30815 0.27536 0.02857
        0.0 0.45977 0.8248
        0.68631 0.76358 0.16174
        0.18869 0.62684 0.95733
        0.18315 0.45826 0.44885
        0.5 0.68637 0.03103
        0.18315 0.46514 0.13804
        0.68274 0.85999 0.75005
        0.31369 0.66583 0.85119
        0.69185 0.4194 0.74694
        0.5 0.11482 0.20554
        0.0 0.10661 0.03195
        0.81685 0.32732 0.11855
        0.0 0.94003 0.71874
        0.5 0.0118 0.57791
        0.19226 0.05841 0.28669
        0.875 0.87263 0.21193
        0.0 0.66756 0.96677
        0.625 0.57781 0.12098
        0.5 0.47797 0.88582
        0.0 0.5541 0.5748
        0.19226 0.94159 0.78669
        0.0 0.73217 0.88476
        0.81131 0.27155 0.31382
        0.0 0.49738 0.9498
        0.625 0.60867 0.78094
        0.30815 0.7081 0.33594
        0.5 0.9882 0.07791
        0.69185 0.58061 0.24694
        0.31369 0.35634 0.56623
        0.0 0.02469 0.03303
        0.125 0.81812 0.39353
        0.19226 0.18462 0.16481
        0.30815 0.4194 0.74694
        0.875 0.18188 0.89353
        0.68631 0.35634 0.56623
        0.81131 0.14239 0.42977
        0.81685 0.46514 0.13804
        0.0 0.29393 0.58751
        0.69185 0.72464 0.52857
        0.19226 0.86248 0.48239
        0.80774 0.05023 0.67578
        0.19226 0.13752 0.98239
        0.5 0.46928 0.01082
        0.18315 0.54175 0.94885
        0.0 0.78979 0.06978
        0.5 0.88518 0.70554
        0.68274 0.93675 0.44262
        0.0 0.26783 0.38476
        0.125 0.26749 0.20522
        0.69185 0.51166 0.0545
        0.68631 0.33418 0.35119
        0.31369 0.64367 0.06623
        0.0 0.80895 0.95314
        0.5 0.52203 0.38582
        0.18315 0.32732 0.11855
        0.5 0.31735 0.00369
        0.31726 0.80967 0.57074
        0.80774 0.99347 0.97598
        0.31369 0.13055 0.52802
        0.5 0.3907 0.07432
        0.30815 0.58061 0.24694
        0.625 0.22024 0.77086
        0.0 0.00409 0.65803
        0.81685 0.45826 0.44885
        0.18869 0.85761 0.92977
        0.19226 0.00653 0.47598
        0.18315 0.67268 0.61855
        0.81131 0.17647 0.62122
        0.0 0.63022 0.32054
        0.19226 0.81538 0.66481
        0.18869 0.14239 0.42977
        0.18315 0.53487 0.63804
        0.5 0.20902 0.38817
        0.0 0.36978 0.82054
        0.81685 0.7311 0.41934
        0.0 0.50263 0.4498
        0.5 0.98214 0.70291
        0.19226 0.94977 0.17578
        0.68631 0.66583 0.85119
        0.80774 0.94159 0.78669
        0.31726 0.19034 0.07074
        0.30815 0.61146 0.6569
        0.68274 0.80967 0.57074
        0.375 0.77977 0.27086
        0.5 0.60931 0.57432
        0.18869 0.17647 0.62122
        0.68274 0.06326 0.94262
        0.0 0.8662 0.64443
        0.625 0.92339 0.33538
        0.5 0.68266 0.50369
        0.31726 0.94933 0.0617
        0.18869 0.82353 0.12122
        0.5 0.79098 0.88817
        0.81131 0.85761 0.92977
        0.68631 0.64367 0.06623
        0.31369 0.81289 0.84136
        0.31726 0.85999 0.75005
        0.875 0.9355 0.89549
        0.5 0.23917 0.59191
        0.19226 0.05023 0.67578
        0.625 0.42219 0.62098
        0.0 0.1338 0.14443
        0.31726 0.14001 0.25005
        0.18315 0.26891 0.91934
        0.625 0.07662 0.83538
        0.0 0.05997 0.21874
        0.31369 0.18712 0.34136
        0.69185 0.48834 0.5545
        0.5 0.93475 0.51059
        0.5 0.7016 0.9021
        0.0 0.70607 0.08751
        0.81685 0.4877 0.76304
        0.5 0.41462 0.82825
        0.31726 0.98708 0.25994
        0.375 0.39134 0.28094
        0.875 0.26749 0.20522
        0.30815 0.47271 0.34339
        0.18315 0.5123 0.26304
        0.19226 0.99347 0.97598
        0.81131 0.72845 0.81382
        0.0 0.97531 0.53303
        0.81685 0.26891 0.91934
        0.81131 0.82353 0.12122
        0.5 0.85271 0.58802
        0.30815 0.38854 0.1569
        0.30815 0.2919 0.83594
        0.69185 0.7081 0.33594
        0.18315 0.6492 0.26202
        0.0 0.54023 0.3248
        0.375 0.92339 0.33538
        0.18869 0.31287 0.6455
        0.31726 0.06326 0.94262
        0.30815 0.52729 0.84339
        0.375 0.60867 0.78094
        0.80774 0.00653 0.47598
        0.5 0.66743 0.38262
        0.68631 0.23643 0.66174
        0.81131 0.31287 0.6455
        0.31369 0.33418 0.35119
        0.31726 0.01292 0.75994
        0.80774 0.81538 0.66481
        0.125 0.73251 0.70522
        0.18869 0.37317 0.45733
        0.68274 0.01292 0.75994
        0.0 0.89339 0.53195
        0.30815 0.72464 0.52857
        0.0 0.4459 0.0748
        0.68631 0.86945 0.02802
        0.80774 0.94977 0.17578
        0.81685 0.6492 0.26202
        0.68631 0.81289 0.84136
        0.68274 0.14001 0.25005
        0.80774 0.86248 0.48239
        0.69185 0.38854 0.1569
        0.375 0.57781 0.12098
        0.69185 0.2919 0.83594
        0.5 0.31363 0.53103
        0.125 0.12737 0.71193
        0.25 0.1522 0.80357
        0.81685 0.5123 0.26304
        0.375 0.07662 0.83538
        0.5 0.06525 0.01059
        0.18315 0.4877 0.76304
        0.5 0.01786 0.20291
        0.31369 0.76358 0.16174
        0.68631 0.18712 0.34136
        0.30815 0.48834 0.5545
        0.81685 0.53487 0.63804
        0.18869 0.68713 0.1455
        0.5 0.76083 0.09191
        0.125 0.0645 0.39549
        0.0 0.67698 0.42673
        0.68274 0.05067 0.5617
        0.25 0.84781 0.30357
        0.0 0.32302 0.92673
        0.125 0.87263 0.21193
        0.5 0.17157 0.51054
        0.5 0.14729 0.08802
        0.75 0.1522 0.80357
        0.81131 0.62684 0.95733
        0.875 0.73251 0.70522
        0.80774 0.05841 0.28669
        0.68274 0.98708 0.25994
        0.69185 0.27536 0.02857
        0.69185 0.61146 0.6569
        0.31726 0.93675 0.44262
        0.81131 0.37317 0.45733
        0.18315 0.3508 0.76202
        0.0 0.33245 0.46677
        0.375 0.22024 0.77086
        0.68274 0.94933 0.0617
        0.68631 0.13055 0.52802
        0.81685 0.54175 0.94885
        0.875 0.81812 0.39353
        0.0 0.19105 0.45314
        0.625 0.77977 0.27086
        0.875 0.12737 0.71193
        0.125 0.18188 0.89353
        0.18315 0.7311 0.41934
        0.80774 0.18462 0.16481
        0.75 0.84781 0.30357
        0.5 0.82844 0.01054
        0.125 0.9355 0.89549
        0.0 0.36153 0.056
        0.81685 0.3508 0.76202
        0.5 0.58539 0.32825
        0.69185 0.52729 0.84339
        0.31726 0.05067 0.5617
        0.5 0.53073 0.51082
        0.80774 0.13752 0.98239
        0.0 0.99591 0.15803
        0.18869 0.27155 0.31382
        0.875 0.0645 0.39549
        0.375 0.42219 0.62098
        0.0 0.63847 0.556
        0.31369 0.23643 0.66174
        0.81131 0.68713 0.1455
        0.69185 0.47271 0.34339
        0.0 0.21022 0.56978
        0.625 0.39134 0.28094
        0.31369 0.86945 0.02802
        0.5 0.29841 0.4021
        0.68274 0.19034 0.07074
        0.30815 0.51166 0.0545
        """

        self.coord = "relative"

        self.cages = """
        16 0.0 -0.39705 0.77884
        12 0.0 -0.16634 0.79524
        15 0.0 -0.05757 1.34533
        12 -0.25475 0.25397 1.48247
        14 0.5 0.81079 0.41895
        12 0.0 0.08992 0.55238
        16 0.5 0.26883 0.19874
        12 0.5 0.48261 0.20147
        14 0.23259 -0.59433 0.94854
        12 -0.5 -0.67027 0.71446
        14 -0.23259 -0.59433 0.94854
        16 0.5 0.05932 0.389
        14 0.23259 0.59433 0.44854
        14 0.5 0.119 0.68822
        12 0.5 0.67027 0.21446
        12 0.0 -0.08992 1.05238
        12 -0.26903 -0.06556 0.6172
        12 0.0 -0.23774 0.54806
        12 -0.26903 0.06556 1.1172
        12 -0.25475 -0.25397 0.98247
        12 0.5 0.41685 0.44471
        12 0.0 0.16634 0.29524
        16 0.0 0.39705 0.27884
        12 0.0 0.23774 0.04806
        14 0.0 -0.243 1.26179
        16 -0.5 -0.05932 0.889
        14 -0.5 -0.119 1.18822
        12 0.25475 0.25397 0.48247
        14 -0.23259 0.59433 1.44854
        14 0.0 0.243 0.76179
        14 -0.5 -0.81079 0.91895
        12 0.25475 -0.25397 0.98247
        12 0.26903 -0.06556 0.6172
        12 -0.5 -0.48261 0.70147
        15 0.0 0.57521 0.12329
        15 0.0 -0.57521 0.62329
        16 -0.5 -0.26883 0.69874
        12 0.26903 0.06556 0.1172
        15 0.0 0.05757 0.84533
        12 -0.5 -0.41685 0.94471
        """

        self.bondlen = 3

        self.cell = """
        14.137381004011147 36.753944148764376 25.74260187174466
        """

        self.density = 0.5094965497744576

        self.cell = cellvectors(a=14.137381004011147,
                                b=36.753944148764376,
                                c=25.74260187174466)
