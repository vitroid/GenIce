# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (18,24,12,0,)
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
        147 181
        234 270
        15 108
        33 259
        220 49
        35 261
        234 57
        297 91
        146 11
        241 107
        61 30
        158 289
        185 210
        211 154
        110 177
        141 268
        135 175
        50 249
        264 58
        70 122
        120 176
        0 235
        110 151
        214 200
        82 92
        291 213
        103 141
        65 154
        216 31
        47 30
        69 77
        249 258
        294 215
        294 213
        166 174
        1 269
        155 178
        290 252
        34 134
        117 96
        178 291
        12 270
        61 55
        240 232
        262 199
        241 233
        263 276
        0 248
        9 82
        17 306
        303 274
        70 124
        181 269
        123 202
        172 100
        224 303
        19 266
        87 159
        27 272
        26 13
        286 3
        218 32
        36 222
        135 28
        190 308
        116 92
        257 47
        95 150
        50 131
        144 73
        304 289
        286 296
        149 53
        190 100
        157 225
        52 185
        173 117
        69 162
        20 237
        1 300
        88 142
        156 154
        131 23
        185 2
        129 304
        241 30
        271 137
        212 244
        134 307
        107 290
        72 199
        42 276
        262 215
        127 260
        198 18
        191 37
        285 208
        249 158
        56 284
        157 253
        228 239
        173 90
        264 31
        119 265
        61 93
        166 90
        203 255
        78 269
        254 18
        46 21
        255 74
        164 108
        165 110
        226 288
        164 202
        163 55
        284 222
        280 82
        303 253
        217 145
        155 222
        197 22
        10 76
        259 104
        64 7
        196 121
        192 208
        153 158
        144 143
        34 35
        14 13
        41 235
        271 37
        52 134
        287 87
        206 170
        140 306
        244 276
        290 47
        122 141
        131 150
        76 55
        166 125
        293 282
        8 85
        197 54
        300 240
        72 57
        32 99
        19 136
        175 207
        194 285
        221 248
        89 195
        239 274
        187 301
        296 246
        266 300
        147 136
        287 75
        9 195
        7 211
        16 298
        253 293
        24 67
        156 74
        75 37
        49 89
        153 58
        300 39
        196 81
        35 180
        236 61
        204 175
        120 208
        120 66
        189 56
        124 126
        33 164
        200 268
        41 221
        40 224
        298 116
        2 180
        137 246
        299 58
        78 279
        210 272
        226 45
        83 274
        27 261
        218 45
        63 135
        80 232
        68 138
        28 23
        157 78
        254 67
        284 112
        68 98
        91 5
        96 169
        229 144
        41 280
        139 12
        152 177
        192 277
        308 67
        176 122
        194 278
        59 98
        281 140
        40 9
        279 250
        68 117
        192 145
        100 273
        111 37
        116 251
        296 159
        281 158
        64 236
        309 288
        86 270
        190 260
        192 260
        42 132
        294 109
        281 267
        272 36
        14 163
        43 38
        186 232
        171 67
        147 121
        25 133
        299 162
        145 176
        302 153
        138 113
        219 309
        219 102
        204 160
        152 191
        115 238
        126 21
        112 91
        206 136
        15 96
        223 279
        53 123
        133 29
        174 138
        218 4
        59 163
        155 102
        309 5
        16 239
        63 88
        118 167
        284 309
        7 200
        310 60
        254 277
        172 79
        52 240
        10 193
        161 31
        0 130
        210 307
        275 97
        212 188
        10 236
        226 168
        228 63
        51 207
        150 311
        212 52
        310 86
        79 177
        216 246
        93 252
        278 71
        20 62
        143 23
        128 233
        81 188
        189 238
        201 12
        285 101
        278 145
        203 301
        205 301
        129 258
        205 65
        248 183
        115 244
        178 91
        100 3
        43 30
        90 13
        179 119
        49 228
        48 227
        64 203
        109 60
        277 101
        196 39
        296 151
        127 278
        2 244
        129 51
        48 116
        137 162
        195 183
        257 187
        194 66
        118 60
        45 5
        297 222
        187 74
        139 53
        134 39
        66 214
        77 289
        25 103
        311 142
        280 293
        210 56
        109 148
        119 44
        164 169
        40 282
        184 228
        29 46
        239 92
        302 283
        71 243
        267 209
        160 227
        97 148
        149 97
        33 8
        124 305
        275 270
        167 237
        43 174
        176 305
        302 69
        205 290
        249 28
        225 136
        43 245
        275 106
        156 124
        151 283
        71 214
        299 289
        184 83
        272 263
        219 57
        24 127
        14 245
        55 38
        63 143
        172 101
        114 188
        217 133
        32 199
        6 251
        257 236
        2 112
        34 80
        186 242
        72 5
        259 54
        165 273
        99 215
        195 6
        10 107
        99 106
        166 54
        198 101
        306 209
        120 29
        144 150
        171 3
        68 76
        46 305
        155 218
        179 197
        48 143
        118 213
        253 248
        121 266
        256 114
        193 128
        217 146
        203 268
        80 256
        85 149
        181 295
        161 153
        17 69
        103 46
        7 187
        184 130
        78 182
        161 258
        50 77
        4 288
        11 214
        256 132
        264 137
        104 169
        271 273
        173 245
        20 97
        154 21
        217 254
        229 142
        25 194
        123 215
        216 283
        247 159
        17 247
        146 243
        259 201
        84 183
        299 159
        189 297
        230 250
        152 260
        227 89
        271 286
        8 44
        34 94
        95 129
        295 250
        105 283
        212 242
        160 23
        308 208
        193 98
        50 306
        229 298
        179 15
        220 298
        135 311
        64 65
        225 224
        223 225
        111 110
        62 106
        109 202
        19 295
        286 287
        24 79
        47 14
        211 11
        193 252
        125 38
        27 115
        140 207
        85 22
        112 168
        22 86
        204 229
        171 191
        302 267
        245 265
        168 115
        227 84
        185 261
        76 26
        65 93
        111 231
        87 105
        56 4
        8 113
        292 62
        291 62
        184 48
        73 88
        0 282
        280 251
        81 206
        202 86
        196 80
        220 130
        89 92
        42 188
        29 11
        170 182
        177 3
        107 255
        103 243
        73 51
        138 265
        83 251
        262 234
        118 199
        171 231
        35 238
        205 126
        94 263
        211 70
        41 223
        269 206
        54 113
        197 169
        224 235
        157 295
        24 285
        287 231
        147 230
        20 294
        308 133
        256 181
        79 231
        216 209
        66 141
        190 198
        151 58
        221 182
        21 268
        191 277
        274 183
        262 139
        275 123
        252 163
        26 174
        161 247
        186 81
        303 82
        266 279
        240 114
        307 276
        104 125
        292 167
        156 200
        139 60
        94 132
        207 142
        128 13
        310 237
        105 162
        127 18
        168 36
        119 104
        96 113
        18 243
        117 38
        247 304
        25 198
        242 94
        233 74
        292 288
        172 75
        42 39
        59 265
        108 53
        152 273
        160 73
        282 6
        180 297
        70 301
        182 230
        146 122
        17 246
        93 233
        36 4
        16 88
        125 59
        238 307
        57 106
        33 148
        261 232
        16 84
        175 131
        235 83
        40 230
        220 6
        165 264
        111 105
        148 12
        72 291
        28 51
        201 22
        84 130
        71 305
        121 170
        19 114
        258 209
        293 250
        180 263
        1 132
        221 9
        165 75
        241 26
        102 99
        90 98
        173 179
        189 45
        44 201
        49 204
        140 304
        102 213
        219 167
        108 44
        223 170
        226 178
        1 186
        281 311
        95 267
        27 242
        126 255
        310 149
        128 257
        234 237
        95 77
        87 31
        15 85
        32 292
        """

        self.waters = """
        9e-05 0.83257 0.84758
        0.0 0.625 0.9576
        0.625 0.0 0.04214
        0.41752 9e-05 0.57448
        0.0 0.5 0.08713
        0.547 0.20381 0.10346
        0.625 0.25 0.83333
        0.0 0.5 0.42046
        0.58257 0.58248 0.24115
        0.37491 0.29244 0.85677
        0.7962 0.5932 0.35906
        0.125 0.625 0.46231
        0.5 0.125 0.20435
        0.875 0.125 0.33359
        0.125 0.25 0.33307
        0.0 0.5 0.2462
        0.375 0.875 0.79564
        0.25 0.125 0.6664
        0.83248 0.99991 0.51424
        0.74991 0.04244 0.93045
        0.625 0.625 0.16666
        0.5462 0.3432 0.43679
        0.7962 0.203 0.22987
        0.5 0.5 0.75379
        0.875 0.375 0.53768
        0.625 0.375 0.5
        0.75 0.875 0.33333
        0.797 0.57881 0.02599
        0.74991 0.70748 0.73622
        0.29244 0.99991 0.47656
        0.32881 0.6568 0.34948
        0.797 0.5932 0.64093
        0.125 0.5 0.12898
        0.45381 0.797 0.22987
        0.25 0.125 0.00026
        0.422 0.20381 0.02625
        0.95748 0.70757 0.06956
        0.078 0.1568 0.58583
        0.20381 0.797 0.3076
        0.2182 0.4212 0.97401
        0.25 0.125 0.87102
        0.375 0.5 0.87102
        0.3432 0.6712 0.98386
        0.32881 0.672 0.31719
        0.41744 0.41752 0.24115
        0.41752 0.41744 0.09219
        0.5 0.375 0.46231
        0.20381 0.4068 0.35906
        0.375 0.5 0.79564
        0.25 0.125 0.79564
        0.375 0.375 0.7088
        0.99991 0.70748 0.73622
        0.875 0.125 0.99974
        0.125 0.5 0.20435
        0.70757 0.95748 0.26378
        0.125 0.875 0.33359
        0.04252 0.29244 0.06956
        0.70748 0.08257 0.14323
        0.578 0.7818 0.64041
        0.4212 0.203 0.30734
        0.16744 0.16752 0.18091
        0.20381 0.7818 0.35958
        0.70748 0.70757 0.14323
        0.5932 0.797 0.77012
        0.04244 0.74991 0.40289
        0.29244 0.99991 0.40289
        0.91744 0.62491 0.47656
        0.125 0.625 0.53768
        0.7962 0.578 0.30708
        0.125 0.25 0.66692
        0.20381 0.4068 0.43679
        0.70757 9e-05 0.47656
        0.5 0.125 0.12898
        0.99991 0.58248 0.75885
        0.70757 9e-05 0.40289
        0.922 0.8432 0.58583
        0.875 0.75 0.33307
        0.2182 0.422 0.69291
        9e-05 0.58257 0.90782
        0.797 0.3432 0.5632
        0.4068 0.20381 0.97427
        0.625 0.625 0.95812
        0.99991 0.16744 0.84758
        0.375 0.625 0.83333
        0.99991 0.70748 0.80989
        0.7962 0.453 0.22987
        0.875 0.125 0.20435
        0.0 0.625 0.62479
        0.3432 0.797 0.77012
        0.37509 0.29252 0.80989
        0.7962 0.203 0.3076
        0.5 0.0 0.08713
        9e-05 0.16752 0.81908
        0.375 0.0 0.37573
        0.125 0.875 0.99974
        0.0 0.375 0.70906
        0.95757 0.70748 0.26378
        0.70757 0.70748 0.1901
        0.6712 0.328 0.31719
        0.08248 0.70757 0.14323
        0.547 0.3432 0.5632
        0.875 0.75 0.53768
        0.125 0.875 0.12898
        0.70757 0.62509 0.47656
        0.29244 0.04252 0.26378
        0.0 0.375 0.62479
        0.83248 0.83257 0.15242
        0.625 0.625 0.37521
        0.20381 0.547 0.22987
        0.29244 0.91752 0.1901
        0.70748 0.99991 0.59711
        0.95748 0.24991 0.59711
        0.70748 0.95757 0.06956
        0.70757 0.70748 0.26378
        0.625 0.0 0.95787
        0.625 0.625 0.04239
        9e-05 0.29252 0.80989
        0.0 0.625 0.2912
        0.16752 0.16744 0.15242
        0.29244 0.29252 0.26378
        0.16744 0.99991 0.48575
        0.24991 0.29244 0.93045
        0.125 0.25 0.46231
        0.08257 0.70748 0.1901
        0.45381 0.6568 0.43679
        0.375 0.0 0.29146
        0.5 0.5 0.42046
        0.70748 0.99991 0.52343
        0.7962 0.2182 0.35958
        0.0 0.625 0.70906
        0.99991 0.83248 0.81908
        0.25009 0.29252 0.73622
        0.2182 0.7962 0.97375
        0.375 0.0 0.5
        0.125 0.25 0.0
        0.8432 0.922 0.74749
        0.8432 0.9212 0.91917
        0.375 0.375 0.62453
        0.625 0.625 0.29093
        0.29244 0.29252 0.1901
        0.375 0.0 0.70854
        0.875 0.75 0.46231
        0.5 0.0 0.75379
        0.3432 0.547 0.77012
        9e-05 0.41752 0.75885
        0.375 0.625 0.5
        0.08257 0.37509 0.47656
        0.1568 0.07881 0.91917
        0.5 0.875 0.20435
        0.875 0.5 0.20435
        9e-05 0.29252 0.73622
        0.625 0.0 0.62427
        0.453 0.6568 0.5632
        0.75 0.875 0.6664
        0.41744 0.99991 0.42552
        0.203 0.7962 0.10346
        0.58257 9e-05 0.42552
        0.3432 0.7962 0.89654
        0.5932 0.797 0.69239
        0.203 0.7818 0.64067
        0.6568 0.453 0.77012
        0.875 0.75 0.66692
        0.203 0.4068 0.64093
        0.25 0.125 0.33333
        0.20381 0.797 0.22987
        0.70748 0.74991 0.59711
        0.625 0.0 0.29146
        0.91752 0.29244 0.14323
        0.70748 0.70757 0.06956
        0.07881 0.922 0.2525
        0.5 0.5 0.91287
        0.203 0.6568 0.5632
        0.797 0.5932 0.5632
        0.0 0.375 0.2912
        0.57881 0.797 0.30734
        0.1568 0.078 0.74749
        0.29244 0.37491 0.47656
        0.58248 0.99991 0.57448
        0.453 0.7962 0.10346
        0.04244 0.29252 0.26378
        0.375 0.0 0.04214
        0.25009 0.95757 0.93045
        0.6568 0.45381 0.89654
        0.0 0.625 0.83333
        0.62491 0.70748 0.80989
        0.797 0.20381 0.02573
        0.7818 0.57881 0.97401
        0.95757 0.25009 0.40289
        0.5932 0.7962 0.97427
        0.29252 0.29244 0.06956
        0.5 0.375 0.53768
        0.203 0.4068 0.5632
        0.29252 0.37509 0.52343
        0.6712 0.3432 0.34948
        0.75 0.375 0.5
        0.625 0.375 0.83333
        0.375 0.375 0.95812
        0.9212 0.078 0.2525
        0.70748 0.62491 0.52343
        0.29252 0.29244 0.14323
        0.7962 0.3432 0.43679
        0.5462 0.203 0.22987
        0.125 0.875 0.20435
        0.9212 0.8432 0.41416
        0.4068 0.203 0.77012
        0.29244 0.24991 0.40289
        0.75009 0.70757 0.93045
        0.24991 0.95748 0.73622
        0.16752 9e-05 0.51424
        0.6568 0.328 0.68281
        0.0 0.375 0.04188
        0.20381 0.6568 0.43679
        0.75 0.875 0.00026
        0.29252 0.91744 0.14323
        0.875 0.375 0.46231
        0.25 0.625 0.16667
        0.672 0.3432 0.65052
        0.25 0.625 0.5
        0.203 0.5462 0.10346
        0.875 0.125 0.12898
        0.37491 0.08248 0.80989
        0.625 0.5 0.87102
        0.078 0.9212 0.08083
        0.3432 0.5462 0.89654
        0.75 0.875 0.87102
        0.5932 0.7962 0.89654
        0.58248 0.58257 0.09219
        0.625 0.5 0.79564
        0.75 0.875 0.79564
        0.6568 0.203 0.77012
        0.4068 0.20381 0.89654
        0.0 0.5 0.57953
        0.6568 0.32881 0.98386
        0.625 0.0 0.37573
        0.625 0.25 0.16666
        0.62509 0.70757 0.85677
        0.0 0.625 0.37547
        0.75 0.375 0.16667
        0.375 0.375 0.04239
        0.62509 0.91752 0.80989
        0.7818 0.20381 0.97375
        0.57881 0.7818 0.35933
        0.875 0.75 0.0
        0.83257 9e-05 0.48575
        0.578 0.7962 0.02625
        0.20381 0.422 0.30708
        0.422 0.2182 0.64041
        0.125 0.875 0.66666
        9e-05 0.70757 0.85677
        0.625 0.625 0.7088
        0.6568 0.20381 0.89654
        0.0 0.375 0.83333
        0.4212 0.2182 0.35933
        0.375 0.875 0.87102
        0.08248 0.37491 0.52343
        0.70757 0.75009 0.40289
        0.375 0.0 0.95787
        0.0 0.375 0.37547
        0.7818 0.578 0.69291
        0.5 0.0 0.2462
        0.5 0.625 0.53768
        0.672 0.32881 0.01615
        0.375 0.375 0.16666
        0.203 0.7962 0.02573
        0.625 0.625 0.62453
        0.375 0.375 0.29093
        0.99991 0.29244 0.93045
        0.7818 0.203 0.69266
        0.7962 0.5932 0.43679
        9e-05 0.70757 0.93045
        0.70757 0.08248 0.1901
        0.29252 0.25009 0.59711
        0.0 0.625 0.04188
        0.5 0.5 0.57953
        0.375 0.75 0.83333
        0.83257 0.83248 0.18091
        0.328 0.6712 0.01615
        0.125 0.25 0.53768
        0.625 0.0 0.5
        0.99991 0.41744 0.90782
        0.99991 0.29244 0.85677
        0.625 0.0 0.70854
        0.37509 0.08257 0.85677
        0.797 0.2182 0.64067
        0.922 0.07881 0.08083
        0.91752 0.62509 0.52343
        0.29252 9e-05 0.59711
        0.04252 0.75009 0.59711
        0.797 0.45381 0.10346
        0.3432 0.672 0.68281
        0.375 0.375 0.37521
        0.5 0.875 0.12898
        0.875 0.5 0.12898
        0.625 0.125 0.87102
        0.375 0.75 0.16666
        0.5 0.0 0.91287
        0.375 0.0 0.62427
        0.29252 0.04244 0.06956
        0.625 0.125 0.79564
        0.328 0.6568 0.65052
        0.0 0.375 0.9576
        0.07881 0.1568 0.41416
        0.875 0.125 0.66666
        0.62491 0.91744 0.85677
        0.2182 0.797 0.69266
        0.5 0.625 0.46231
        0.4068 0.203 0.69239
        0.203 0.4212 0.02599
        0.29252 9e-05 0.52343
        0.797 0.20381 0.10346
        0.91744 0.29252 0.1901
        0.75009 0.04252 0.73622
        """

        self.coord = "relative"

        self.cages = """
        14 0.17009 0.16974 0.10968
        14 0.83026 0.00035 0.44301
        14 0.82991 0.83026 0.10968
        14 0.16974 0.17009 0.22365
        14 0.99965 0.82991 0.77634
        14 0.16974 0.99965 0.44301
        14 0.00035 0.83026 0.89032
        14 0.83026 0.82991 0.22365
        14 0.82991 0.99965 0.55698
        14 0.00035 0.17009 0.77634
        14 0.99965 0.16974 0.89032
        14 0.17009 0.00035 0.55698
        12 0.5 0.0 0.0
        12 0.0 0.5 0.33333
        12 0.5 0.5 0.66666
        15 0.5 0.0 0.46324
        15 0.0 0.5 0.79657
        15 0.0 0.5 0.87009
        15 0.5 0.5 0.1299
        15 0.5 0.0 0.53676
        15 0.5 0.5 0.20342
        12 0.0 0.0 0.16667
        12 0.0 0.0 0.5
        12 0.0 0.0 0.83333
        12 0.5 0.0 0.33437
        12 0.0 0.5 0.6677
        12 0.0 0.5 0.99896
        12 0.5 0.5 0.00103
        12 0.5 0.0 0.66563
        12 0.5 0.5 0.33229
        12 0.5 0.0 0.16666
        12 0.0 0.5 0.49999
        12 0.0 0.5 0.16667
        12 0.5 0.5 0.83332
        12 0.5 0.0 0.83334
        12 0.5 0.5 0.5
        15 0.0 0.0 0.03938
        15 0.0 0.0 0.37271
        15 0.0 0.0 0.29395
        15 0.0 0.0 0.70604
        15 0.0 0.0 0.96062
        15 0.0 0.0 0.62728
        14 0.31199 0.68478 0.06458
        14 0.31522 0.62721 0.39791
        14 0.68801 0.31522 0.06458
        14 0.68478 0.31199 0.26875
        14 0.37279 0.68801 0.73124
        14 0.68478 0.37279 0.39791
        14 0.62721 0.31522 0.93542
        14 0.31522 0.68801 0.26875
        14 0.68801 0.37279 0.60208
        14 0.62721 0.31199 0.73124
        14 0.37279 0.68478 0.93542
        14 0.31199 0.62721 0.60208
        """

        self.bondlen = 3

        self.cell = """
        22.08 0.0 0.0
        -11.039999999999994 19.121840915560405 0.0
        1.5024824996270077e-14 2.6023760268370632e-14 245.37401325395933
        """

        self.density = 0.09001792793959411

        self.cell = cellvectors(a=22.08,
                                b=22.08,
                                c=245.37401325395933,
                                C=119.99999999999999)
