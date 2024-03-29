# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (32,16,16,4,)
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
        166 216
        125 280
        17 364
        25 374
        158 223
        46 269
        156 62
        335 254
        234 120
        352 20
        45 129
        149 201
        263 136
        121 212
        256 255
        68 258
        309 192
        22 130
        149 241
        59 63
        48 122
        41 138
        232 307
        274 71
        209 232
        234 351
        214 324
        172 359
        367 89
        353 188
        231 107
        195 35
        40 203
        279 163
        227 81
        231 322
        19 118
        210 189
        200 140
        229 273
        182 237
        118 302
        193 187
        161 164
        308 292
        92 352
        364 80
        7 224
        361 102
        293 282
        236 24
        246 359
        50 262
        310 113
        203 152
        197 323
        318 71
        25 141
        207 230
        364 365
        288 185
        0 258
        301 340
        158 345
        260 141
        346 35
        31 120
        329 224
        272 31
        174 57
        274 280
        35 365
        210 254
        210 253
        98 199
        90 89
        249 357
        305 102
        277 137
        6 349
        200 111
        204 233
        32 374
        24 130
        131 343
        286 380
        114 99
        179 53
        143 320
        24 13
        162 201
        329 319
        330 319
        34 228
        262 104
        386 215
        132 379
        259 158
        42 383
        214 191
        2 363
        16 51
        83 20
        167 53
        133 307
        337 47
        166 276
        301 163
        19 144
        106 264
        143 382
        204 64
        146 280
        265 378
        256 165
        347 246
        104 359
        78 302
        287 14
        30 158
        68 227
        279 26
        301 82
        10 62
        179 216
        119 350
        239 104
        125 251
        343 328
        293 370
        319 154
        133 279
        240 151
        28 154
        293 378
        218 105
        121 271
        288 244
        269 223
        51 337
        324 275
        32 385
        56 302
        20 70
        1 3
        211 74
        197 51
        2 106
        285 220
        369 15
        339 173
        76 384
        74 144
        221 39
        70 328
        229 226
        204 237
        127 181
        33 205
        252 159
        372 15
        14 328
        186 198
        80 268
        153 85
        300 245
        284 257
        66 243
        325 244
        0 182
        126 55
        189 268
        101 69
        157 154
        149 383
        147 55
        214 100
        142 135
        29 228
        74 213
        124 283
        306 84
        298 215
        314 242
        362 90
        110 383
        137 21
        285 350
        124 140
        291 123
        300 11
        234 93
        28 267
        77 202
        324 5
        173 196
        268 190
        203 65
        117 350
        285 83
        344 104
        129 307
        40 96
        325 49
        50 327
        49 176
        178 295
        29 194
        160 348
        148 192
        311 59
        16 331
        147 372
        45 209
        314 143
        0 331
        181 62
        2 331
        37 375
        25 371
        329 107
        116 141
        148 261
        1 26
        298 198
        166 376
        86 254
        85 67
        129 250
        150 145
        266 86
        248 9
        165 5
        101 320
        337 345
        56 3
        127 79
        368 177
        259 321
        42 255
        236 63
        239 76
        149 79
        40 68
        240 52
        4 244
        34 113
        222 191
        186 95
        235 276
        86 321
        28 224
        218 159
        252 379
        358 99
        184 188
        132 115
        247 26
        10 200
        286 213
        240 144
        241 375
        123 178
        164 112
        338 9
        211 57
        12 350
        29 208
        16 147
        217 110
        129 342
        360 57
        225 61
        226 326
        352 128
        262 347
        161 207
        92 353
        60 111
        53 376
        25 136
        123 9
        94 371
        118 279
        209 185
        304 49
        80 321
        386 349
        225 307
        44 248
        236 73
        131 171
        90 43
        91 78
        200 333
        126 345
        303 379
        61 122
        54 250
        151 213
        209 39
        16 243
        330 11
        8 267
        362 169
        18 243
        361 356
        113 112
        270 96
        257 172
        13 218
        369 103
        151 366
        153 217
        278 176
        381 98
        54 84
        190 280
        18 64
        14 64
        157 136
        118 52
        32 10
        284 302
        289 251
        266 369
        135 107
        370 58
        264 242
        114 213
        132 22
        256 62
        165 383
        6 132
        150 277
        18 296
        284 122
        355 268
        55 66
        100 85
        57 117
        30 335
        155 290
        91 52
        385 102
        229 346
        151 72
        183 205
        114 332
        184 299
        116 333
        330 136
        147 103
        180 217
        48 340
        13 67
        1 347
        27 272
        241 192
        92 299
        340 26
        245 145
        46 313
        178 358
        33 341
        1 52
        250 175
        139 22
        180 255
        224 102
        59 349
        7 322
        108 216
        97 75
        148 127
        70 353
        8 263
        314 269
        315 75
        103 339
        167 215
        325 84
        88 194
        125 355
        12 168
        19 174
        338 299
        304 179
        190 365
        305 322
        121 95
        291 315
        157 333
        275 261
        142 341
        133 50
        298 159
        60 267
        338 336
        17 86
        309 191
        354 194
        30 196
        306 36
        342 84
        54 221
        370 135
        81 323
        114 199
        354 71
        355 69
        39 175
        36 221
        296 81
        64 315
        359 174
        109 365
        88 230
        121 342
        182 47
        288 175
        42 130
        356 154
        45 271
        269 321
        380 220
        286 171
        387 295
        201 140
        46 355
        195 138
        94 351
        297 310
        288 306
        182 227
        99 171
        14 381
        309 115
        92 171
        127 255
        339 331
        309 24
        256 31
        126 196
        142 378
        311 73
        335 373
        17 326
        265 169
        124 351
        146 109
        19 327
        125 138
        275 31
        160 378
        305 38
        296 363
        34 362
        146 228
        150 90
        344 247
        370 322
        249 166
        174 56
        53 349
        124 219
        211 128
        225 163
        137 134
        293 297
        155 220
        266 313
        88 43
        168 128
        301 39
        79 100
        137 43
        287 97
        27 5
        197 242
        373 173
        44 199
        208 35
        260 219
        214 241
        60 157
        289 382
        133 344
        265 89
        4 276
        153 139
        0 65
        142 367
        282 134
        27 111
        61 306
        131 220
        148 120
        319 245
        67 252
        80 253
        78 172
        233 152
        178 75
        46 15
        329 367
        272 10
        291 152
        273 195
        199 20
        170 110
        205 88
        206 43
        177 117
        360 3
        48 36
        71 238
        87 207
        155 74
        66 363
        265 134
        312 253
        314 254
        328 233
        177 290
        326 316
        323 196
        361 107
        366 290
        179 379
        222 85
        295 188
        185 281
        168 72
        77 281
        7 282
        8 283
        377 228
        372 223
        184 233
        311 376
        41 69
        318 377
        287 358
        168 360
        226 109
        300 58
        101 15
        203 237
        294 78
        33 206
        126 96
        153 63
        229 164
        42 192
        23 204
        65 243
        303 59
        49 198
        263 219
        240 360
        131 248
        291 381
        6 191
        354 183
        334 202
        308 202
        63 170
        325 235
        99 9
        95 281
        289 316
        187 238
        235 95
        144 177
        318 273
        119 336
        231 141
        343 299
        79 217
        263 317
        225 175
        312 101
        267 385
        120 375
        323 47
        363 258
        312 259
        41 316
        37 201
        135 145
        156 162
        261 115
        334 311
        239 257
        312 335
        156 374
        346 238
        285 211
        258 75
        367 206
        304 159
        7 11
        115 67
        348 341
        190 238
        245 89
        338 332
        165 162
        111 162
        352 332
        292 281
        73 215
        382 189
        305 260
        270 237
        327 368
        226 187
        119 290
        8 371
        176 54
        30 242
        167 357
        41 253
        143 313
        40 66
        5 375
        160 87
        34 348
        300 134
        212 185
        364 138
        32 116
        72 119
        248 381
        231 330
        326 210
        343 83
        387 152
        377 164
        334 105
        344 384
        366 91
        368 91
        77 198
        11 317
        212 76
        303 105
        38 58
        187 189
        50 82
        297 33
        6 170
        354 87
        222 386
        4 250
        108 176
        382 69
        83 336
        261 100
        292 276
        123 184
        93 156
        60 374
        87 169
        297 277
        324 110
        13 139
        21 169
        161 205
        2 197
        150 341
        294 3
        38 317
        342 384
        146 193
        208 230
        275 180
        93 371
        186 278
        21 230
        247 56
        251 313
        362 194
        183 348
        98 387
        105 216
        207 346
        183 377
        347 368
        208 112
        339 47
        82 232
        17 251
        4 108
        38 356
        181 140
        170 130
        273 316
        55 106
        232 76
        366 380
        294 246
        287 353
        29 274
        51 223
        271 278
        351 272
        327 172
        160 310
        310 161
        28 317
        73 249
        373 266
        277 58
        186 167
        12 294
        356 145
        386 252
        333 219
        303 139
        292 45
        283 37
        94 385
        103 345
        48 239
        12 380
        23 188
        358 98
        23 97
        373 320
        247 122
        308 278
        271 221
        218 249
        296 97
        340 262
        337 270
        298 334
        260 94
        282 206
        264 320
        195 274
        27 283
        289 193
        65 315
        108 202
        68 295
        318 193
        72 332
        212 36
        37 93
        77 244
        234 181
        372 264
        96 81
        44 336
        304 357
        82 257
        246 117
        284 163
        61 384
        235 357
        308 376
        18 270
        155 44
        113 21
        70 387
        286 128
        222 236
        116 361
        369 259
        173 106
        109 112
        23 227
        22 180
        """

        self.waters = """
        0.0 0.125 0.2672
        0.125 0.25 0.4828
        0.1875 0.92832 0.23698
        0.19669 0.4375 0.46978
        0.31713 0.94213 0.61484
        0.75 0.25 0.80139
        0.625 0.25 0.72396
        0.5 0.43287 0.93525
        0.69213 0.43287 0.86484
        0.0 0.05787 0.35936
        0.18287 0.07169 0.84113
        0.69669 0.55787 0.92718
        0.0 0.56713 0.43525
        0.3125 0.75 0.71102
        0.625 0.75 0.32838
        0.69213 0.05332 0.17718
        0.8125 0.92832 0.23698
        0.0 0.375 0.1258
        0.625 0.75 0.27604
        0.5 0.875 0.47396
        0.5 0.375 0.37421
        0.82169 0.5625 0.01302
        0.3125 0.07169 0.73698
        0.125 0.57169 0.30131
        0.5 0.75 0.72396
        0.5 0.75 0.875
        0.19669 0.25 0.50876
        0.875 0.25 0.82838
        0.69669 0.25 0.90419
        0.82169 0.125 0.05131
        0.375 0.57169 0.19869
        0.25 0.25 0.80139
        0.30787 0.06713 0.86484
        0.30332 0.25 0.99124
        0.0 0.875 0.02604
        0.69669 0.44213 0.07283
        0.0 0.43287 0.56476
        0.69213 0.55332 0.82283
        0.0 0.375 0.92162
        0.0 0.06713 0.56476
        0.5 0.25 0.27604
        0.5 0.25 0.125
        0.5 0.875 0.7672
        0.69669 0.25 0.99124
        0.68287 0.05787 0.38516
        0.67832 0.06713 0.59113
        0.81713 0.92832 0.15887
        0.0 0.44669 0.24124
        0.0 0.4375 0.53898
        0.19213 0.43287 0.63516
        0.67832 0.0625 0.51302
        0.875 0.75 0.22396
        0.19669 0.0625 0.46978
        0.625 0.25 0.67162
        0.19669 0.25 0.59582
        0.5 0.125 0.2328
        0.32169 0.5625 0.48698
        0.5 0.5 0.44861
        0.0 0.5 0.94861
        0.875 0.07169 0.69869
        0.69213 0.06713 0.86484
        0.32169 0.625 0.55131
        0.30787 0.94669 0.82283
        0.8125 0.94669 0.71978
        0.75 0.75 0.30139
        0.8125 0.05332 0.28022
        0.5 0.05332 0.25876
        0.1875 0.55332 0.71978
        0.3125 0.25 0.28898
        0.5 0.05787 0.14064
        0.5 0.44669 0.34582
        0.5 0.93287 0.06476
        0.125 0.25 0.41442
        0.68287 0.75 0.68525
        0.5 0.875 0.42162
        0.125 0.92832 0.30131
        0.67832 0.625 0.55131
        0.0 0.75 0.625
        0.0 0.75 0.47396
        0.0 0.75 0.77604
        0.69213 0.56713 0.13516
        0.3125 0.57169 0.26302
        0.80332 0.9375 0.53022
        0.68287 0.44213 0.38516
        0.32169 0.43287 0.59113
        0.0 0.625 0.7328
        0.0 0.44669 0.15419
        0.5 0.75 0.02604
        0.625 0.25 0.0172
        0.69669 0.9375 0.96978
        0.82169 0.0625 0.98698
        0.0 0.9375 0.46102
        0.19213 0.56713 0.36484
        0.5 0.625 0.83558
        0.30787 0.43287 0.86484
        0.68287 0.55787 0.61484
        0.5 0.44669 0.25876
        0.25 0.75 0.30139
        0.5 0.05332 0.34582
        0.19213 0.93287 0.36484
        0.0 0.55332 0.75876
        0.5 0.125 0.16442
        0.30332 0.25 0.90419
        0.8125 0.25 0.21102
        0.67832 0.4375 0.51302
        0.19213 0.94669 0.67718
        0.3125 0.05332 0.21978
        0.30332 0.94213 0.92718
        0.19213 0.06713 0.63516
        0.0 0.625 0.07838
        0.8125 0.07169 0.76302
        0.81713 0.07169 0.84113
        0.0 0.5 0.05139
        0.0 0.625 0.02604
        0.31713 0.05787 0.38516
        0.3125 0.42832 0.73698
        0.18287 0.94213 0.88516
        0.67832 0.375 0.44869
        0.32169 0.9375 0.48698
        0.875 0.25 0.41442
        0.375 0.42832 0.80131
        0.67832 0.43287 0.59113
        0.19669 0.5625 0.53022
        0.0 0.125 0.33558
        0.0 0.55332 0.84582
        0.81713 0.05787 0.11484
        0.5 0.375 0.2328
        0.1875 0.75 0.78898
        0.32169 0.56713 0.40887
        0.5 0.125 0.57838
        0.5 0.94669 0.74124
        0.875 0.75 0.3758
        0.375 0.25 0.72396
        0.5 0.125 0.52604
        0.69669 0.5625 0.96978
        0.17832 0.875 0.94869
        0.69213 0.75 0.89064
        0.82169 0.4375 0.98698
        0.69213 0.25 0.10936
        0.1875 0.94669 0.71978
        0.0 0.75 0.82759
        0.30787 0.75 0.89064
        0.30332 0.9375 0.96978
        0.18287 0.92832 0.15887
        0.5 0.0 0.44861
        0.0 0.0 0.94861
        0.0 0.875 0.07838
        0.6875 0.05332 0.21978
        0.3125 0.55332 0.78022
        0.8125 0.75 0.78898
        0.0 0.125 0.97396
        0.19669 0.05787 0.42718
        0.68287 0.25 0.31476
        0.0 0.875 0.7328
        0.82169 0.06713 0.90887
        0.67832 0.93287 0.40887
        0.5 0.875 0.83558
        0.81713 0.94213 0.88516
        0.625 0.57169 0.19869
        0.19213 0.55332 0.67718
        0.30332 0.75 0.00876
        0.30332 0.4375 0.03022
        0.69213 0.94669 0.82283
        0.19669 0.9375 0.53022
        0.17832 0.375 0.05131
        0.625 0.07169 0.80131
        0.5 0.94669 0.65419
        0.68287 0.42832 0.65887
        0.19669 0.44213 0.42718
        0.69669 0.75 0.00876
        0.6875 0.07169 0.73698
        0.125 0.75 0.3758
        0.80332 0.75 0.49124
        0.1875 0.25 0.21102
        0.5 0.625 0.47396
        0.19669 0.94213 0.57283
        0.125 0.25 0.62421
        0.67832 0.125 0.44869
        0.19213 0.05332 0.32283
        0.375 0.25 0.67162
        0.1875 0.07169 0.76302
        0.18287 0.75 0.81476
        0.0 0.375 0.2672
        0.30332 0.0625 0.03022
        0.0 0.375 0.33558
        0.875 0.75 0.58558
        0.80787 0.43287 0.63516
        0.30332 0.75 0.09582
        0.19213 0.44669 0.32283
        0.375 0.75 0.12421
        0.69669 0.75 0.09582
        0.6875 0.42832 0.73698
        0.5 0.625 0.7672
        0.17832 0.93287 0.09113
        0.69669 0.0625 0.03022
        0.625 0.25 0.08558
        0.3125 0.44669 0.21978
        0.125 0.75 0.22396
        0.0 0.55787 0.64064
        0.5 0.125 0.37421
        0.0 0.94669 0.84582
        0.81713 0.75 0.81476
        0.0 0.94213 0.64064
        0.6875 0.25 0.28898
        0.875 0.57169 0.30131
        0.375 0.25 0.0172
        0.5 0.25 0.97396
        0.5 0.5625 0.03898
        0.82169 0.375 0.05131
        0.80332 0.94213 0.57283
        0.30787 0.56713 0.13516
        0.5 0.625 0.42162
        0.80332 0.55787 0.57283
        0.32169 0.93287 0.40887
        0.8125 0.42832 0.76302
        0.80787 0.55332 0.67718
        0.31713 0.07169 0.65887
        0.0 0.94669 0.75876
        0.31713 0.75 0.68525
        0.0 0.625 0.87421
        0.80332 0.75 0.40419
        0.0 0.25 0.57759
        0.8125 0.55332 0.71978
        0.75 0.75 0.19861
        0.5 0.25 0.92242
        0.32169 0.875 0.55131
        0.17832 0.56713 0.09113
        0.1875 0.44669 0.28022
        0.0 0.0 0.05139
        0.30332 0.44213 0.07283
        0.69669 0.4375 0.03022
        0.375 0.75 0.91442
        0.67832 0.875 0.55131
        0.80787 0.44669 0.32283
        0.30787 0.55332 0.82283
        0.5 0.625 0.6258
        0.6875 0.75 0.71102
        0.8125 0.44669 0.28022
        0.5 0.75 0.07759
        0.80332 0.5625 0.53022
        0.32169 0.125 0.44869
        0.6875 0.55332 0.78022
        0.25 0.75 0.19861
        0.6875 0.92832 0.26302
        0.19213 0.75 0.60936
        0.82169 0.875 0.94869
        0.80332 0.4375 0.46978
        0.32169 0.4375 0.51302
        0.80787 0.93287 0.36484
        0.5 0.75 0.67242
        0.32169 0.06713 0.59113
        0.0 0.125 0.1258
        0.125 0.42832 0.69869
        0.5 0.44213 0.14064
        0.18287 0.57169 0.15887
        0.3125 0.94669 0.78022
        0.375 0.07169 0.80131
        0.875 0.75 0.5172
        0.1875 0.05332 0.28022
        0.69213 0.44669 0.17718
        0.18287 0.55787 0.88516
        0.1875 0.42832 0.76302
        0.80332 0.25 0.50876
        0.81713 0.55787 0.88516
        0.375 0.92832 0.19869
        0.625 0.75 0.9828
        0.0 0.25 0.17242
        0.625 0.25 0.8758
        0.625 0.75 0.12421
        0.875 0.75 0.17162
        0.6875 0.57169 0.26302
        0.80332 0.25 0.59582
        0.125 0.25 0.82838
        0.375 0.25 0.08558
        0.69669 0.05787 0.07283
        0.125 0.25 0.77604
        0.5 0.875 0.6258
        0.0 0.375 0.97396
        0.875 0.25 0.62421
        0.32169 0.0625 0.51302
        0.82169 0.93287 0.09113
        0.80787 0.75 0.60936
        0.5 0.4375 0.96102
        0.81713 0.42832 0.84113
        0.125 0.75 0.5172
        0.67832 0.56713 0.40887
        0.19669 0.75 0.40419
        0.375 0.75 0.32838
        0.125 0.75 0.58558
        0.18287 0.05787 0.11484
        0.80332 0.05787 0.42718
        0.80787 0.05332 0.32283
        0.68287 0.94213 0.61484
        0.30332 0.5625 0.96978
        0.0 0.5625 0.46102
        0.31713 0.25 0.31476
        0.375 0.75 0.27604
        0.17832 0.4375 0.98698
        0.0 0.625 0.66442
        0.0 0.44213 0.35936
        0.82169 0.625 0.94869
        0.0 0.0625 0.53898
        0.19669 0.75 0.49124
        0.125 0.07169 0.69869
        0.31713 0.42832 0.65887
        0.17832 0.43287 0.90887
        0.19669 0.55787 0.57283
        0.5 0.0 0.55139
        0.80787 0.06713 0.63516
        0.5 0.55332 0.74124
        0.17832 0.5625 0.01302
        0.80787 0.94669 0.67718
        0.5 0.375 0.16442
        0.0 0.05332 0.15419
        0.125 0.75 0.17162
        0.875 0.92832 0.30131
        0.30787 0.25 0.10936
        0.82169 0.43287 0.90887
        0.30332 0.05787 0.07283
        0.69669 0.94213 0.92718
        0.30787 0.05332 0.17718
        0.81713 0.57169 0.15887
        0.30332 0.55787 0.92718
        0.1875 0.57169 0.23698
        0.875 0.25 0.77604
        0.31713 0.55787 0.61484
        0.18287 0.44213 0.11484
        0.67832 0.9375 0.48698
        0.68287 0.57169 0.34113
        0.5 0.06713 0.93525
        0.625 0.75 0.91442
        0.0 0.05332 0.24124
        0.19213 0.25 0.39064
        0.0 0.875 0.87421
        0.0 0.875 0.66442
        0.30787 0.44669 0.17718
        0.80787 0.25 0.39064
        0.8125 0.57169 0.23698
        0.0 0.25 0.375
        0.0 0.25 0.22396
        0.0 0.25 0.52604
        0.17832 0.0625 0.98698
        0.5 0.375 0.57838
        0.80787 0.56713 0.36484
        0.5 0.375 0.52604
        0.6875 0.44669 0.21978
        0.5 0.56713 0.06476
        0.875 0.25 0.4828
        0.17832 0.9375 0.01302
        0.75 0.25 0.69861
        0.80332 0.44213 0.42718
        0.18287 0.42832 0.84113
        0.31713 0.44213 0.38516
        0.31713 0.57169 0.34113
        0.5 0.9375 0.03898
        0.69213 0.93287 0.13516
        0.0 0.125 0.92162
        0.5 0.55332 0.65419
        0.31713 0.92832 0.34113
        0.67832 0.5625 0.48698
        0.32169 0.375 0.44869
        0.17832 0.06713 0.90887
        0.82169 0.9375 0.01302
        0.3125 0.92832 0.26302
        0.81713 0.44213 0.11484
        0.82169 0.56713 0.09113
        0.0 0.93287 0.43525
        0.5 0.0625 0.96102
        0.80332 0.0625 0.46978
        0.81713 0.25 0.18525
        0.17832 0.625 0.94869
        0.5 0.55787 0.85936
        0.625 0.92832 0.19869
        0.18287 0.25 0.18525
        0.5 0.94213 0.85936
        0.625 0.42832 0.80131
        0.68287 0.07169 0.65887
        0.17832 0.125 0.05131
        0.375 0.75 0.9828
        0.25 0.25 0.69861
        0.0 0.75 0.42242
        0.68287 0.92832 0.34113
        0.30787 0.93287 0.13516
        0.6875 0.94669 0.78022
        0.5 0.5 0.55139
        0.375 0.25 0.8758
        0.875 0.42832 0.69869
        0.5 0.25 0.32759
        """

        self.coord = "relative"

        self.cages = """
        12 0.0 0.98148 0.39691
        12 -0.25 0.25 -0.75
        15 0.0 -0.75 -0.26703
        15 0.0 0.75 0.33436
        12 0.0 -0.98148 -0.39691
        12 0.5 0.75 -0.05177
        15 0.5 1.25 0.83436
        12 0.73148 0.75 0.64691
        12 1.0 1.5 1.5
        16 0.5 -0.75 0.625
        12 0.23148 0.25 0.14691
        12 0.5 1.0 1.0
        12 0.0 0.25 0.30177
        12 0.0 -0.75 -0.55177
        12 0.25 -0.25 0.75
        12 0.26852 0.25 0.35309
        14 0.5 -0.53674 0.69792
        15 0.0 -0.75 -0.33436
        12 0.0 1.0 0.5
        12 -0.26852 0.25 -0.64691
        16 0.0 -0.75 -0.125
        14 1.0 1.03674 0.80208
        14 0.78674 0.75 0.05208
        15 0.5 0.25 0.41564
        12 -0.23148 -0.25 -0.14691
        14 0.78674 1.25 0.94792
        14 0.0 -0.53674 -0.19792
        12 0.26852 -0.25 0.64691
        16 -0.5 0.75 0.375
        14 0.0 0.53674 0.19792
        12 1.0 1.51852 1.39691
        12 0.5 1.25 1.05177
        14 0.28674 -0.75 -0.44792
        15 0.0 0.75 0.26703
        12 -0.25 -0.25 -0.25
        12 0.5 -0.25 0.80177
        15 0.0 -0.25 -0.01703
        14 -0.78674 -1.25 -0.94792
        12 0.5 1.01852 0.89691
        15 0.0 -0.25 -0.08436
        15 0.5 0.25 0.48297
        14 1.0 -0.03674 1.19792
        12 0.5 0.98148 0.10309
        14 0.5 -0.03674 0.30208
        12 0.0 -0.25 -0.30177
        15 0.5 0.75 0.23297
        14 -0.28674 0.75 0.44792
        14 -0.28674 1.25 0.55208
        15 0.5 -0.25 0.51703
        15 0.0 0.25 0.08436
        14 0.5 1.03674 0.69792
        12 0.5 0.25 0.19823
        12 0.25 0.25 0.25
        12 0.0 0.48148 -0.39691
        15 0.5 0.75 0.16564
        12 0.5 0.51852 0.10309
        14 1.28674 1.75 1.44792
        15 0.5 1.25 0.76703
        12 0.5 1.48148 0.89691
        12 0.5 1.5 1.0
        12 0.23148 -0.25 -0.14691
        14 0.21326 1.25 0.94792
        14 0.5 0.53674 0.30208
        12 0.76852 0.25 1.14691
        15 0.5 -0.25 0.58436
        12 0.0 0.75 0.55177
        15 0.0 0.25 0.01703
        16 0.0 0.75 0.125
        """

        self.bondlen = 3

        self.cell = """
        13.350439626923315 13.350439626923315 109.77461386228596
        """

        self.density = 0.5927484501192363

        self.cell = cellvectors(a=13.350439626923315,
                                b=13.350439626923315,
                                c=109.77461386228596)
