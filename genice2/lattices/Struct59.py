# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (24,16,16,0,)
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
        270 208
        108 157
        129 252
        49 192
        318 219
        199 78
        262 73
        143 201
        182 232
        135 89
        32 123
        219 228
        93 298
        55 167
        208 13
        286 212
        296 104
        161 35
        40 250
        124 185
        228 245
        135 283
        292 37
        57 312
        289 98
        5 297
        74 288
        92 90
        275 281
        282 277
        150 4
        264 245
        188 110
        26 240
        293 172
        223 158
        220 195
        233 265
        203 245
        312 155
        150 33
        260 206
        86 24
        263 122
        25 91
        47 166
        100 98
        255 252
        79 70
        107 38
        36 228
        59 236
        309 34
        121 62
        289 282
        29 75
        228 112
        59 117
        105 61
        268 245
        115 158
        319 8
        41 111
        267 142
        25 56
        300 53
        193 206
        100 69
        311 127
        235 56
        278 171
        215 42
        137 6
        42 138
        272 226
        284 96
        170 309
        255 234
        217 259
        73 95
        84 276
        273 303
        243 202
        310 66
        12 43
        238 303
        27 185
        283 203
        301 222
        199 212
        270 150
        105 267
        70 239
        18 103
        249 48
        51 248
        77 165
        241 13
        28 112
        214 221
        169 60
        292 162
        70 125
        153 280
        15 187
        16 187
        261 63
        139 304
        99 44
        142 165
        78 282
        272 236
        75 173
        181 79
        29 211
        72 151
        204 30
        264 318
        175 248
        146 169
        115 78
        269 285
        51 251
        284 93
        237 123
        52 124
        269 153
        108 175
        90 69
        171 185
        16 225
        260 112
        197 195
        178 265
        111 17
        175 302
        77 95
        294 293
        233 232
        88 90
        72 81
        319 162
        37 126
        261 107
        50 117
        291 186
        254 197
        128 304
        241 147
        301 189
        261 313
        52 299
        56 244
        54 246
        94 145
        76 2
        46 159
        102 88
        219 82
        268 38
        244 247
        268 235
        120 237
        209 20
        249 231
        303 298
        169 39
        2 187
        57 117
        186 6
        157 151
        263 216
        148 67
        22 41
        8 173
        142 103
        36 110
        257 202
        2 138
        137 98
        176 24
        289 186
        131 82
        269 75
        116 195
        47 280
        189 183
        174 37
        129 225
        242 154
        136 313
        164 84
        227 89
        152 20
        80 318
        242 313
        15 129
        119 0
        293 226
        254 315
        121 133
        175 62
        199 74
        215 63
        119 51
        23 190
        96 101
        10 234
        136 118
        239 49
        12 62
        287 114
        106 258
        31 101
        159 206
        63 224
        48 290
        105 205
        300 280
        160 216
        279 239
        244 220
        30 101
        265 191
        253 84
        311 158
        95 307
        213 183
        305 224
        3 295
        121 31
        210 23
        131 277
        214 85
        73 310
        42 220
        229 71
        114 198
        97 309
        182 200
        40 235
        243 80
        141 117
        189 19
        174 109
        40 220
        28 193
        10 210
        11 209
        7 240
        9 65
        272 138
        0 180
        92 240
        279 98
        278 315
        268 259
        1 104
        115 186
        163 198
        26 68
        57 14
        10 180
        274 64
        91 242
        316 14
        79 274
        27 315
        25 265
        96 229
        267 181
        274 61
        9 126
        46 191
        28 110
        223 317
        281 125
        94 21
        271 55
        178 67
        76 124
        46 218
        314 311
        11 21
        33 42
        233 206
        307 30
        81 216
        247 305
        205 132
        249 181
        207 134
        26 106
        192 68
        100 277
        79 192
        278 299
        301 104
        154 207
        279 297
        194 34
        43 184
        37 248
        212 276
        306 147
        229 298
        10 71
        144 178
        270 259
        88 192
        302 43
        21 287
        29 184
        255 122
        146 166
        208 266
        243 314
        249 165
        119 101
        210 122
        262 253
        99 5
        237 222
        222 156
        145 202
        120 86
        137 49
        300 253
        139 290
        156 89
        111 39
        85 65
        72 298
        41 97
        63 244
        68 246
        77 285
        193 44
        285 231
        118 23
        25 13
        91 247
        238 184
        38 208
        190 251
        177 85
        152 64
        144 200
        36 306
        270 89
        189 59
        284 285
        271 147
        176 148
        238 252
        138 230
        316 160
        168 109
        130 123
        271 112
        159 145
        149 230
        197 247
        163 188
        277 74
        151 134
        204 92
        154 305
        66 231
        3 127
        67 232
        143 139
        3 290
        319 303
        316 50
        209 132
        210 134
        271 82
        106 54
        107 56
        199 164
        135 217
        300 310
        260 38
        80 94
        284 133
        204 54
        263 213
        94 24
        33 250
        309 308
        131 6
        11 198
        125 64
        9 22
        194 179
        209 125
        120 218
        29 146
        27 215
        257 281
        72 190
        61 103
        128 87
        314 286
        73 173
        295 280
        19 256
        168 62
        130 200
        16 174
        6 170
        16 315
        227 13
        230 195
        103 304
        232 196
        251 171
        293 116
        314 140
        273 81
        143 84
        205 297
        291 113
        119 71
        299 207
        143 295
        257 99
        214 39
        7 69
        17 288
        261 266
        255 155
        76 302
        136 316
        126 211
        127 212
        134 216
        318 113
        27 2
        294 230
        170 74
        188 275
        28 5
        83 166
        105 88
        66 18
        67 19
        264 24
        219 282
        184 162
        154 160
        182 313
        60 126
        241 191
        107 233
        257 140
        40 156
        144 242
        194 295
        301 148
        58 287
        136 224
        238 45
        292 187
        22 54
        64 140
        100 102
        218 178
        47 34
        260 241
        59 296
        55 275
        202 114
        222 172
        108 251
        70 205
        14 23
        1 217
        46 114
        76 149
        83 177
        55 5
        172 183
        224 185
        22 83
        129 294
        85 240
        137 308
        90 246
        253 111
        161 252
        108 124
        273 161
        294 35
        53 18
        215 4
        286 317
        118 171
        144 256
        0 168
        12 211
        49 167
        4 266
        130 196
        91 266
        152 311
        272 15
        275 239
        256 213
        149 225
        4 197
        53 276
        250 104
        278 109
        141 213
        11 44
        221 7
        53 3
        267 92
        279 110
        153 201
        269 310
        57 122
        9 39
        190 71
        139 165
        132 87
        264 217
        128 20
        227 203
        262 201
        221 34
        250 259
        281 198
        236 35
        174 12
        273 312
        97 288
        32 141
        258 133
        8 229
        116 156
        243 58
        83 211
        283 113
        15 155
        258 231
        58 113
        226 296
        306 203
        0 133
        80 317
        297 102
        45 75
        254 52
        226 33
        304 276
        145 44
        35 183
        93 180
        214 166
        52 305
        235 227
        146 153
        118 207
        120 148
        78 317
        223 20
        223 21
        152 290
        163 191
        160 256
        302 292
        291 82
        287 86
        308 7
        31 65
        188 147
        182 50
        41 47
        306 291
        36 289
        128 286
        179 170
        263 161
        262 169
        299 157
        135 237
        115 179
        204 65
        164 17
        151 180
        159 196
        179 127
        32 296
        51 109
        319 155
        58 158
        1 123
        194 164
        258 30
        193 163
        141 200
        254 149
        77 201
        168 157
        176 196
        99 132
        102 167
        1 176
        308 68
        106 177
        96 95
        45 234
        66 307
        32 172
        131 167
        274 48
        87 140
        43 225
        288 69
        150 116
        130 218
        221 17
        60 162
        121 177
        31 248
        61 87
        48 18
        81 14
        50 19
        283 86
        26 181
        312 236
        60 173
        97 246
        45 93
        142 307
        8 234
        """

        self.waters = """
        0.375 0.28896 0.25
        0.68603 0.67761 0.64944
        0.81103 0.44138 0.90588
        0.31103 0.04018 0.60116
        0.5 0.56143 0.99296
        0.68897 0.87334 0.19102
        0.31103 0.94138 0.90588
        0.0 0.06143 0.99296
        0.68603 0.30022 0.55898
        0.68897 0.20637 0.93398
        0.5 0.34233 0.40955
        0.0 0.84233 0.40955
        0.31397 0.30022 0.94102
        0.31397 0.69978 0.05898
        0.81103 0.45982 0.39884
        0.625 0.4243 0.75
        0.5 0.39078 0.9094
        0.81397 0.07228 0.83634
        0.5 0.08561 0.44454
        0.0 0.58561 0.44454
        0.0 0.93857 0.49296
        0.0 0.83036 0.53455
        0.5 0.17294 0.95138
        0.68603 0.42772 0.33634
        0.875 0.75659 0.625
        0.18897 0.65981 0.1424
        0.18603 0.10269 0.15662
        0.68603 0.46929 0.00352
        0.81397 0.82239 0.14944
        0.25 0.25988 0.75
        0.625 0.18995 0.25
        0.875 0.24342 0.125
        0.5 0.5898 0.60468
        0.68603 0.57228 0.83634
        0.18603 0.07228 0.83634
        0.125 0.48529 0.66014
        0.0 0.83036 0.96546
        0.68603 0.30022 0.94102
        0.68603 0.69978 0.05898
        0.81397 0.17761 0.85057
        0.0 0.62528 0.86796
        0.5 0.12528 0.86796
        0.81103 0.54018 0.89884
        0.18897 0.34019 0.85761
        0.81397 0.82239 0.35057
        0.31397 0.30022 0.55898
        0.375 0.74342 0.375
        0.31103 0.12666 0.80898
        0.31103 0.05862 0.40588
        0.31103 0.95982 0.10116
        0.81103 0.55862 0.40588
        0.68603 0.32239 0.14944
        0.18897 0.45982 0.10116
        0.5 0.07463 0.56954
        0.5 0.14646 0.07638
        0.5 0.87472 0.13204
        0.0 0.64646 0.07638
        0.68603 0.46929 0.49648
        0.31103 0.84019 0.6424
        0.81103 0.54018 0.60116
        0.75 0.25988 0.75
        0.625 0.01471 0.33986
        0.18897 0.29364 0.06602
        0.81103 0.55862 0.09412
        0.31103 0.95982 0.39884
        0.81397 0.19978 0.05898
        0.5 0.14646 0.42363
        0.0 0.64646 0.42363
        0.31103 0.05862 0.09412
        0.81397 0.03072 0.99648
        0.125 0.96217 0.25
        0.68603 0.32239 0.35057
        0.0 0.37472 0.36796
        0.68897 0.20637 0.56602
        0.625 0.98529 0.83986
        0.375 0.25659 0.625
        0.0 0.4144 0.94454
        0.0 0.16964 0.46546
        0.875 0.9243 0.75
        0.25 0.01844 0.25
        0.68897 0.84019 0.6424
        0.0 0.42537 0.43046
        0.5 0.85354 0.92363
        0.31103 0.20637 0.93398
        0.81397 0.07228 0.66366
        0.0 0.16964 0.03455
        0.125 0.75659 0.625
        0.68897 0.95982 0.39884
        0.625 0.01471 0.16014
        0.31397 0.67761 0.85057
        0.68897 0.05862 0.09412
        0.31397 0.60269 0.15662
        0.81397 0.10269 0.15662
        0.18897 0.29364 0.43398
        0.81397 0.80022 0.55898
        0.81397 0.19978 0.44102
        0.875 0.24342 0.375
        0.5 0.07463 0.93046
        0.0 0.93857 0.00704
        0.68897 0.87334 0.30898
        0.81397 0.96929 0.00352
        0.75 0.24013 0.25
        0.68897 0.95982 0.10116
        0.68897 0.05862 0.40588
        0.81103 0.62666 0.69102
        0.75 0.01844 0.25
        0.31103 0.15981 0.1424
        0.81103 0.65981 0.1424
        0.0 0.37472 0.13204
        0.5 0.34233 0.09046
        0.0 0.84233 0.09046
        0.68897 0.12666 0.80898
        0.68897 0.79364 0.06602
        0.375 0.81005 0.75
        0.31103 0.79364 0.43398
        0.125 0.9243 0.75
        0.31397 0.57228 0.83634
        0.68603 0.53072 0.50352
        0.625 0.46217 0.25
        0.625 0.28896 0.25
        0.18897 0.70637 0.56602
        0.125 0.24342 0.125
        0.5 0.43857 0.49296
        0.5 0.65767 0.59046
        0.0 0.42537 0.06954
        0.18603 0.92772 0.33634
        0.625 0.25659 0.875
        0.375 0.98529 0.66014
        0.81397 0.96929 0.49648
        0.375 0.4243 0.75
        0.5 0.66964 0.46546
        0.5 0.9144 0.94454
        0.81397 0.92772 0.33634
        0.25 0.24013 0.25
        0.31397 0.42772 0.33634
        0.375 0.71104 0.75
        0.75 0.51844 0.25
        0.18603 0.96929 0.00352
        0.875 0.48529 0.83986
        0.0 0.06143 0.50704
        0.5 0.92537 0.43046
        0.5 0.56143 0.50704
        0.81397 0.10269 0.34338
        0.0 0.0898 0.60468
        0.31397 0.60269 0.34338
        0.68897 0.79364 0.43398
        0.125 0.21104 0.75
        0.31103 0.79364 0.06602
        0.0 0.67294 0.54863
        0.18897 0.44138 0.90588
        0.5 0.5898 0.89532
        0.18897 0.37334 0.30898
        0.18603 0.96929 0.49648
        0.18603 0.17761 0.64944
        0.25 0.51844 0.25
        0.68603 0.39732 0.65662
        0.18897 0.62666 0.80898
        0.18897 0.37334 0.19102
        0.18603 0.89732 0.65662
        0.625 0.74342 0.375
        0.125 0.51471 0.33986
        0.18897 0.44138 0.59412
        0.875 0.31005 0.75
        0.125 0.78896 0.25
        0.875 0.03784 0.75
        0.0 0.10923 0.4094
        0.18603 0.17761 0.85057
        0.5 0.92537 0.06954
        0.31397 0.32239 0.14944
        0.875 0.21104 0.75
        0.375 0.98529 0.83986
        0.68603 0.42772 0.16366
        0.31397 0.57228 0.66366
        0.625 0.25659 0.625
        0.5 0.33036 0.96546
        0.0 0.32707 0.04863
        0.81103 0.70637 0.56602
        0.18603 0.19978 0.05898
        0.18897 0.65981 0.35761
        0.25 0.98156 0.75
        0.31397 0.32239 0.35057
        0.125 0.0757 0.25
        0.68603 0.60269 0.34338
        0.18897 0.54018 0.60116
        0.125 0.31005 0.75
        0.81103 0.45982 0.10116
        0.18603 0.89732 0.84338
        0.68603 0.39732 0.84338
        0.18603 0.82239 0.14944
        0.0 0.57463 0.56954
        0.81103 0.37334 0.30898
        0.25 0.74013 0.25
        0.375 0.01471 0.16014
        0.875 0.78896 0.25
        0.125 0.03784 0.75
        0.18897 0.54018 0.89884
        0.68603 0.69978 0.44102
        0.31397 0.53072 0.99648
        0.18603 0.82239 0.35057
        0.75 0.98156 0.75
        0.5 0.60923 0.4094
        0.0 0.15767 0.59046
        0.5 0.82707 0.45138
        0.125 0.75659 0.875
        0.68897 0.15981 0.1424
        0.875 0.96217 0.25
        0.75 0.74013 0.25
        0.375 0.46217 0.25
        0.5 0.66964 0.03455
        0.0 0.9102 0.39532
        0.5 0.4102 0.39532
        0.375 0.25659 0.875
        0.625 0.98529 0.66014
        0.31397 0.53072 0.50352
        0.0 0.15767 0.90955
        0.68603 0.53072 0.99648
        0.18897 0.45982 0.39884
        0.625 0.71104 0.75
        0.31397 0.69978 0.44102
        0.68897 0.84019 0.85761
        0.0 0.57463 0.93046
        0.0 0.0898 0.89532
        0.18897 0.62666 0.69102
        0.0 0.89078 0.5906
        0.875 0.51471 0.16014
        0.31397 0.39732 0.84338
        0.625 0.53784 0.75
        0.18897 0.70637 0.93398
        0.81397 0.80022 0.94102
        0.81103 0.29364 0.43398
        0.125 0.48529 0.83986
        0.31103 0.15981 0.35761
        0.81103 0.65981 0.35761
        0.875 0.68995 0.25
        0.5 0.33036 0.53455
        0.0 0.67294 0.95138
        0.875 0.48529 0.66014
        0.31397 0.67761 0.64944
        0.18897 0.34019 0.6424
        0.18603 0.92772 0.16366
        0.0 0.10923 0.0906
        0.375 0.74342 0.125
        0.375 0.5757 0.25
        0.5 0.85354 0.57638
        0.0 0.58561 0.05546
        0.875 0.75659 0.875
        0.5 0.08561 0.05546
        0.18897 0.55862 0.09412
        0.81103 0.29364 0.06602
        0.18603 0.10269 0.34338
        0.81103 0.62666 0.80898
        0.81103 0.37334 0.19102
        0.31397 0.39732 0.65662
        0.68897 0.12666 0.69102
        0.31397 0.46929 0.00352
        0.5 0.39078 0.5906
        0.18897 0.55862 0.40588
        0.5 0.87472 0.36796
        0.375 0.18995 0.25
        0.68603 0.67761 0.85057
        0.625 0.74342 0.125
        0.68603 0.60269 0.15662
        0.81397 0.17761 0.64944
        0.31397 0.46929 0.49648
        0.75 0.75988 0.75
        0.125 0.68995 0.25
        0.5 0.60923 0.0906
        0.875 0.0757 0.25
        0.81103 0.70637 0.93398
        0.31103 0.20637 0.56602
        0.5 0.65767 0.90955
        0.5 0.82707 0.04863
        0.75 0.48156 0.75
        0.0 0.4144 0.55546
        0.375 0.01471 0.33986
        0.31103 0.87334 0.19102
        0.68897 0.04018 0.60116
        0.68897 0.94138 0.90588
        0.5 0.4102 0.10468
        0.0 0.9102 0.10468
        0.31103 0.12666 0.69102
        0.31103 0.87334 0.30898
        0.81397 0.89732 0.84338
        0.25 0.75988 0.75
        0.125 0.24342 0.375
        0.18603 0.19978 0.44102
        0.68897 0.94138 0.59412
        0.18603 0.80022 0.55898
        0.68897 0.04018 0.89884
        0.0 0.89078 0.9094
        0.18603 0.03072 0.50352
        0.31103 0.84019 0.85761
        0.81103 0.34019 0.85761
        0.375 0.53784 0.75
        0.25 0.48156 0.75
        0.18603 0.07228 0.66366
        0.68603 0.57228 0.66366
        0.81397 0.92772 0.16366
        0.0 0.32707 0.45138
        0.31397 0.42772 0.16366
        0.5 0.12528 0.63204
        0.0 0.62528 0.63204
        0.0 0.35354 0.92363
        0.0 0.35354 0.57638
        0.81397 0.03072 0.50352
        0.125 0.51471 0.16014
        0.18603 0.80022 0.94102
        0.68897 0.15981 0.35761
        0.18603 0.03072 0.99648
        0.31103 0.04018 0.89884
        0.5 0.17294 0.54863
        0.31103 0.94138 0.59412
        0.81103 0.44138 0.59412
        0.625 0.5757 0.25
        0.5 0.9144 0.55546
        0.5 0.43857 0.00704
        0.875 0.51471 0.33986
        0.81397 0.89732 0.65662
        0.625 0.81005 0.75
        0.81103 0.34019 0.6424
        """

        self.coord = "relative"

        self.cages = """
        12 0.0 0.5 0.0
        14 0.0 0.01492 0.35943
        15 0.5 0.24342 0.08366
        12 0.24411 0.12286 0.51408
        14 0.5 0.48508 -0.14057
        15 0.5 0.67626 0.25
        14 0.5 0.94721 0.25
        12 0.5 0.0 0.5
        15 1.0 0.74342 0.08366
        12 -0.25589 0.37714 0.48592
        12 0.25589 -0.37714 -0.01408
        12 0.0 0.12839 0.75
        12 0.25589 0.37714 0.01408
        14 0.0 -0.2974 0.75
        14 0.5 0.51492 0.35943
        15 -0.5 -0.0967 0.75
        12 0.75589 0.87714 -0.51408
        12 -0.25589 -0.37714 -0.01408
        12 0.5 0.37161 0.25
        12 -0.25589 -0.37714 0.51408
        14 0.5 0.7974 0.25
        15 1.0 1.17626 0.25
        12 0.0 -0.5 0.5
        14 0.0 -0.44721 0.75
        15 0.0 -0.17626 0.75
        12 0.24411 0.12286 -0.01408
        15 1.0 0.5967 0.25
        14 0.0 0.01492 0.14057
        15 0.0 0.25658 0.58366
        12 1.0 0.87161 0.25
        15 -0.5 0.24342 0.41634
        14 0.0 -0.01492 -0.14057
        12 0.24411 0.87714 0.48592
        15 0.5 0.0967 0.25
        15 0.0 0.25658 -0.08366
        14 0.0 -0.01492 0.64057
        12 0.75589 0.12286 -0.01408
        12 -0.5 -0.37161 0.75
        12 0.24411 0.87714 0.01408
        14 0.5 0.51492 0.14057
        12 0.5 1.0 0.0
        14 0.5 0.48508 0.64057
        15 -0.5 -0.24342 -0.08366
        12 -0.25589 0.37714 0.01408
        12 0.75589 0.12286 -0.48592
        15 -0.5 -0.24342 0.58366
        15 0.0 0.74342 0.41634
        12 0.75589 0.87714 0.01408
        14 0.0 0.2974 0.25
        15 0.0 0.4033 0.75
        14 0.0 0.44721 0.25
        14 0.5 0.05279 0.75
        15 -0.5 -0.67626 0.75
        14 0.5 0.2026 0.75
        12 0.25589 -0.37714 -0.48592
        12 0.25589 0.37714 -0.51408
        """

        self.bondlen = 3

        self.cell = """
        13.325063821889895 49.74878164013856 24.403574339604887
        """

        self.density = 0.5912573383843822

        self.cell = cellvectors(a=13.325063821889895,
                                b=49.74878164013856,
                                c=24.403574339604887)
