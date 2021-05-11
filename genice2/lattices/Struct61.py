# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (32,20,8,8,)
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
        366 121
        146 63
        109 336
        370 89
        48 198
        193 102
        294 200
        164 362
        49 341
        289 262
        343 314
        374 55
        296 299
        205 174
        331 34
        362 145
        82 178
        373 119
        77 182
        156 332
        297 17
        9 13
        222 242
        315 379
        139 265
        33 367
        135 13
        24 259
        346 179
        278 246
        161 200
        128 331
        311 121
        147 366
        291 149
        151 125
        215 326
        384 357
        60 69
        238 99
        131 98
        60 65
        122 202
        14 142
        146 5
        5 208
        138 186
        91 93
        170 218
        3 291
        277 130
        250 283
        298 142
        176 126
        334 387
        302 42
        257 331
        247 353
        161 192
        196 59
        197 61
        28 166
        262 7
        166 293
        59 5
        354 164
        27 30
        292 266
        270 268
        10 195
        329 282
        259 138
        318 81
        96 41
        49 170
        50 170
        182 114
        185 305
        94 0
        267 214
        345 134
        308 243
        272 18
        374 92
        358 372
        281 204
        172 328
        258 149
        282 348
        117 293
        331 314
        311 109
        135 81
        253 99
        70 20
        340 384
        173 306
        201 0
        59 58
        321 322
        301 116
        235 350
        229 2
        82 67
        155 381
        9 82
        214 299
        225 187
        233 81
        379 306
        227 224
        346 242
        330 67
        84 246
        359 273
        148 377
        100 324
        206 196
        236 375
        17 380
        83 319
        319 200
        4 17
        206 41
        142 109
        285 359
        99 246
        4 230
        117 105
        239 74
        46 132
        252 34
        14 167
        237 195
        277 270
        130 335
        279 31
        226 169
        345 207
        329 189
        108 287
        281 45
        215 373
        271 261
        79 89
        27 50
        193 111
        105 337
        248 143
        72 281
        15 367
        62 78
        26 65
        15 366
        194 179
        295 119
        225 12
        243 93
        250 343
        34 221
        145 176
        168 344
        284 156
        118 138
        385 199
        237 274
        137 256
        380 128
        116 273
        303 56
        98 305
        286 113
        75 91
        112 288
        20 233
        22 120
        221 376
        135 178
        284 326
        285 328
        295 66
        209 57
        37 265
        257 268
        227 137
        234 127
        208 190
        292 345
        103 90
        275 295
        139 40
        33 142
        67 134
        165 269
        189 372
        116 110
        62 158
        64 25
        342 181
        201 29
        222 129
        104 114
        127 114
        144 182
        351 294
        117 370
        205 195
        162 232
        65 106
        80 365
        103 79
        180 166
        237 153
        43 88
        71 343
        18 189
        4 7
        154 252
        366 151
        173 300
        60 334
        206 204
        218 305
        75 256
        108 110
        290 8
        38 175
        323 212
        169 361
        256 372
        298 327
        378 126
        72 158
        111 319
        43 224
        51 26
        307 207
        300 20
        237 323
        379 376
        369 263
        21 120
        383 13
        174 240
        186 121
        56 307
        290 130
        219 278
        32 58
        85 327
        133 67
        247 276
        349 157
        303 39
        274 69
        226 179
        5 13
        304 23
        88 371
        359 151
        378 115
        317 173
        355 92
        333 183
        81 66
        326 356
        150 239
        266 86
        231 177
        103 57
        33 121
        254 210
        22 353
        65 39
        264 220
        144 298
        340 213
        3 280
        217 158
        8 140
        317 200
        303 6
        225 28
        195 160
        12 211
        220 140
        172 199
        227 56
        236 313
        12 89
        363 48
        203 129
        58 110
        80 39
        243 274
        309 333
        216 371
        253 295
        152 46
        277 212
        338 209
        197 381
        23 94
        211 288
        293 98
        337 263
        108 163
        26 299
        101 313
        101 30
        258 294
        357 191
        373 350
        351 291
        180 177
        235 326
        339 382
        227 165
        22 171
        172 29
        62 191
        38 102
        19 361
        292 378
        383 184
        244 223
        288 310
        131 382
        30 92
        310 100
        38 31
        53 208
        54 209
        112 68
        136 83
        348 213
        15 78
        44 218
        97 106
        11 347
        215 42
        35 153
        132 327
        192 173
        136 232
        24 169
        325 369
        69 160
        21 63
        270 314
        189 80
        323 97
        245 306
        37 52
        87 95
        15 23
        271 318
        219 76
        386 174
        107 233
        47 291
        133 328
        289 164
        311 368
        348 241
        87 335
        144 385
        151 320
        193 47
        154 175
        124 298
        79 280
        234 77
        36 272
        118 357
        203 157
        222 97
        89 112
        363 220
        27 341
        9 255
        95 107
        350 177
        254 70
        286 317
        205 171
        255 7
        334 75
        354 228
        216 90
        16 325
        137 76
        299 207
        258 363
        87 184
        266 6
        159 94
        221 63
        207 261
        101 88
        19 48
        32 208
        61 374
        44 238
        43 288
        90 370
        249 100
        242 96
        14 199
        386 74
        313 47
        304 86
        222 230
        386 302
        84 240
        79 149
        316 132
        176 278
        351 169
        338 148
        271 178
        64 55
        52 336
        301 297
        94 249
        211 101
        36 204
        316 194
        73 318
        186 125
        379 143
        226 3
        235 321
        199 263
        280 377
        329 36
        198 264
        375 312
        10 284
        185 1
        205 381
        84 91
        202 255
        312 302
        84 74
        308 241
        217 320
        354 134
        11 107
        212 267
        219 261
        274 106
        152 361
        141 252
        249 86
        302 25
        83 85
        206 158
        308 244
        16 344
        320 78
        374 90
        22 386
        322 356
        215 239
        155 140
        371 358
        160 240
        141 171
        253 364
        71 248
        345 251
        198 279
        265 123
        297 241
        250 198
        369 168
        360 183
        226 96
        322 55
        258 8
        6 115
        330 287
        347 1
        273 287
        383 7
        105 341
        308 153
        113 114
        120 130
        359 46
        187 356
        128 190
        11 71
        144 46
        364 1
        28 370
        325 77
        168 6
        51 267
        244 21
        340 194
        251 126
        72 357
        244 157
        147 285
        330 289
        27 37
        163 260
        122 242
        286 283
        51 203
        93 157
        210 233
        254 177
        367 159
        301 213
        384 132
        263 115
        346 213
        309 231
        146 223
        375 8
        111 351
        52 368
        126 334
        153 257
        103 148
        325 385
        32 163
        124 360
        296 60
        297 58
        188 88
        324 358
        253 218
        118 327
        16 368
        320 202
        167 100
        50 25
        148 18
        127 162
        17 257
        339 248
        333 245
        183 229
        234 368
        323 268
        184 230
        129 41
        9 78
        349 74
        43 324
        37 211
        354 73
        45 272
        372 387
        204 348
        134 271
        235 332
        266 307
        375 187
        104 283
        342 41
        365 282
        318 119
        35 171
        364 352
        45 338
        16 147
        329 269
        95 300
        77 293
        352 143
        108 29
        155 356
        203 165
        66 245
        122 4
        190 221
        76 307
        154 376
        95 343
        284 214
        115 387
        152 113
        315 2
        133 304
        141 247
        71 128
        25 92
        117 369
        383 107
        275 170
        276 355
        259 111
        2 85
        264 252
        174 64
        340 281
        285 385
        180 321
        382 185
        165 342
        136 317
        313 355
        154 279
        72 259
        124 131
        229 306
        346 365
        10 212
        238 332
        338 194
        122 301
        335 192
        303 68
        216 105
        275 40
        44 321
        316 54
        38 319
        39 387
        377 80
        51 97
        143 190
        225 265
        224 181
        217 96
        85 31
        214 239
        49 98
        196 241
        347 164
        137 324
        245 20
        234 123
        268 230
        264 314
        150 246
        197 220
        52 310
        19 286
        273 202
        365 106
        152 125
        113 232
        363 57
        296 251
        210 73
        352 53
        87 270
        14 159
        131 182
        180 305
        49 360
        10 155
        349 267
        382 104
        304 147
        70 104
        201 260
        146 184
        159 29
        66 352
        249 76
        112 168
        33 118
        236 12
        231 40
        315 279
        315 339
        289 163
        276 61
        11 262
        232 186
        217 24
        355 312
        337 358
        44 64
        311 162
        362 364
        316 361
        290 161
        110 191
        54 48
        223 129
        262 380
        349 120
        201 176
        0 178
        53 260
        292 133
        247 175
        59 62
        83 138
        196 223
        3 181
        19 294
        61 57
        188 272
        35 34
        377 179
        102 209
        269 256
        250 339
        275 333
        290 353
        373 40
        139 50
        175 161
        251 228
        86 344
        367 191
        228 145
        73 350
        35 21
        1 248
        381 55
        68 224
        36 342
        166 123
        0 219
        136 229
        243 282
        188 181
        216 30
        341 336
        160 156
        141 197
        172 378
        371 18
        135 53
        330 255
        162 183
        123 231
        23 82
        332 228
        278 75
        69 91
        296 156
        193 45
        28 322
        344 310
        167 336
        337 167
        70 309
        347 210
        277 140
        145 99
        287 328
        353 312
        236 149
        300 283
        150 119
        276 102
        280 68
        150 261
        124 2
        42 187
        63 335
        260 362
        139 42
        32 380
        254 185
        384 116
        54 31
        24 125
        188 47
        238 240
        56 26
        269 93
        376 192
        127 309
        360 109
        """

        self.waters = """
        0.00813 0.64209 0.77781
        0.10432 0.3515 0.44417
        0.20382 0.13951 0.62156
        0.55741 0.89841 0.06917
        0.35209 0.59808 0.24333
        0.28164 0.5604 0.81573
        0.83545 0.81177 0.22219
        0.26378 0.57304 0.18083
        0.5604 0.21836 0.06573
        0.19308 0.64477 0.99677
        0.64791 0.40192 0.24333
        0.20278 0.44246 0.27802
        0.73623 0.07304 0.06917
        0.21836 0.5604 0.93427
        0.02218 0.86841 0.62177
        0.14377 0.78265 0.90302
        0.9976 0.92778 0.10042
        0.35945 0.56568 0.40281
        0.6485 0.89568 0.55583
        0.37427 0.0926 0.19094
        0.16455 0.31177 0.02781
        0.49177 0.45864 0.7599
        0.59808 0.35209 0.75667
        0.07304 0.73623 0.93083
        0.35209 0.90192 0.00667
        0.77327 0.22472 0.6875
        0.70382 0.63951 0.12156
        0.85209 0.09808 0.74333
        0.80692 0.14477 0.25323
        0.0604 0.71836 0.56573
        0.76378 0.07304 0.68083
        0.3515 0.10432 0.55583
        0.23623 0.57304 0.56917
        0.14791 0.90192 0.74333
        0.41982 0.38669 0.56594
        0.4976 0.42778 0.60042
        0.55755 0.79723 0.72198
        0.85945 0.06568 0.90281
        0.41982 0.11332 0.68406
        0.70937 0.75946 0.3125
        0.92696 0.26378 0.93083
        0.48605 0.73623 0.94417
        0.77528 0.27327 0.9375
        0.75 0.88972 0.875
        0.85523 0.30692 0.50323
        0.50241 0.92778 0.64958
        0.21836 0.9396 0.31573
        0.57223 0.9976 0.85042
        0.4074 0.12573 0.30906
        0.97782 0.13159 0.62177
        0.85523 0.19308 0.74677
        0.62573 0.5926 0.05906
        0.92778 0.9976 0.89958
        0.14209 0.50813 0.72219
        0.39841 0.05741 0.43083
        0.72472 0.22673 0.4375
        0.74055 0.70937 0.0625
        0.56086 0.08541 0.39958
        0.28265 0.64377 0.59698
        0.30692 0.64477 0.75323
        0.75 0.61028 0.375
        0.59808 0.14791 0.49333
        0.27528 0.72673 0.8125
        0.39841 0.44259 0.81917
        0.77327 0.27528 0.5625
        0.68823 0.66455 0.27781
        0.0926 0.37427 0.80906
        0.06568 0.64055 0.15281
        0.70382 0.86049 0.12844
        0.70278 0.55755 0.47198
        0.11028 0.25 0.125
        0.25 0.38972 0.375
        0.3636 0.8636 0.75
        0.00241 0.42778 0.14958
        0.71836 0.4396 0.81573
        0.79618 0.63951 0.62844
        0.86049 0.70382 0.87156
        0.00823 0.04136 0.2599
        0.22472 0.72673 0.9375
        0.64209 0.00813 0.22219
        0.63951 0.79618 0.37156
        0.10159 0.44259 0.93083
        0.09808 0.64791 0.99333
        0.29723 0.05755 0.77802
        0.76378 0.48605 0.69417
        0.27609 0.07119 0.61833
        0.94246 0.79723 0.02802
        0.37427 0.4074 0.05906
        0.70278 0.94246 0.77802
        0.73623 0.01395 0.19417
        0.71836 0.0604 0.43427
        0.72391 0.57119 0.63167
        0.71736 0.14377 0.65302
        0.6485 0.60432 0.69417
        0.01395 0.73623 0.80583
        0.29618 0.36049 0.12156
        0.4396 0.78164 0.06573
        0.55741 0.60159 0.18083
        0.99187 0.14209 0.47219
        0.91459 0.43914 0.60042
        0.89568 0.8515 0.80583
        0.71836 0.02659 0.80927
        0.4976 0.07223 0.64958
        0.63159 0.02218 0.37177
        0.16455 0.18823 0.22219
        0.85791 0.00813 0.52781
        0.60432 0.6485 0.30583
        0.22391 0.42881 0.11833
        0.14477 0.69308 0.50323
        0.04136 0.00823 0.7401
        0.22673 0.72472 0.5625
        0.42778 0.00241 0.85042
        0.77609 0.92881 0.13167
        0.22391 0.07119 0.13167
        0.1485 0.10432 0.19417
        0.86049 0.79618 0.37844
        0.27528 0.77327 0.4375
        0.86841 0.02218 0.37823
        0.23623 0.92696 0.68083
        0.93914 0.41459 0.89958
        0.56086 0.41459 0.85042
        0.14055 0.93433 0.90281
        0.35523 0.69308 0.24677
        0.93433 0.14055 0.09719
        0.12573 0.0926 0.55906
        0.26378 0.92696 0.06917
        0.89568 0.6485 0.44417
        0.08018 0.11332 0.06594
        0.29723 0.44246 0.47198
        0.49187 0.64209 0.97219
        0.50813 0.35791 0.97219
        0.0926 0.12573 0.44094
        0.28164 0.9396 0.43427
        0.02659 0.71836 0.19073
        0.9976 0.57223 0.14958
        0.13159 0.52218 0.87177
        0.25 0.11028 0.875
        0.79063 0.74055 0.8125
        0.28164 0.97342 0.80927
        0.85624 0.21736 0.90302
        0.5604 0.28164 0.18427
        0.52659 0.28164 0.55927
        0.08541 0.93914 0.64958
        0.20382 0.36049 0.62844
        0.13159 0.97782 0.37823
        0.95864 0.50823 0.5099
        0.36841 0.52218 0.87823
        0.06568 0.85945 0.09719
        0.60159 0.94259 0.43083
        0.58541 0.06086 0.10042
        0.86841 0.47782 0.87177
        0.21736 0.85624 0.09698
        0.26378 0.98605 0.19417
        0.5 0.5 0.5
        0.38972 0.25 0.625
        0.64477 0.30692 0.24677
        0.78164 0.47342 0.30927
        0.60159 0.55741 0.81917
        0.35624 0.78265 0.84698
        0.0604 0.78164 0.68427
        0.71836 0.47342 0.44073
        0.42881 0.22391 0.88167
        0.11332 0.08018 0.93406
        0.14791 0.59808 0.50667
        0.07223 0.4976 0.35042
        0.63951 0.70382 0.87844
        0.90192 0.14791 0.25667
        0.94259 0.89841 0.68083
        0.8515 0.89568 0.19417
        0.41459 0.93914 0.10042
        0.9396 0.21836 0.68427
        0.56568 0.35945 0.59719
        0.01395 0.76378 0.44417
        0.29063 0.25946 0.9375
        0.71736 0.35624 0.59698
        0.44246 0.20278 0.72198
        0.94259 0.60159 0.56917
        0.97342 0.28164 0.19073
        0.06086 0.58541 0.89958
        0.49187 0.85791 0.27781
        0.92696 0.23623 0.31917
        0.60432 0.8515 0.94417
        0.10159 0.05741 0.31917
        0.10432 0.1485 0.80583
        0.35791 0.50813 0.02781
        0.07119 0.27609 0.38167
        0.21836 0.97342 0.94073
        0.72673 0.22472 0.0625
        0.61332 0.91982 0.81594
        0.66455 0.81177 0.52781
        0.27609 0.42881 0.63167
        0.22673 0.77528 0.6875
        0.36049 0.29618 0.87844
        0.5 0.0 0.75
        0.43914 0.91459 0.39958
        0.64055 0.43433 0.40281
        0.40192 0.64791 0.75667
        0.57304 0.23623 0.43083
        0.36049 0.20382 0.37156
        0.00813 0.85791 0.47219
        0.39568 0.1485 0.94417
        0.02218 0.63159 0.62823
        0.27327 0.72472 0.1875
        0.5926 0.62573 0.94094
        0.47342 0.78164 0.69073
        0.6364 0.3636 0.5
        0.42696 0.73623 0.81917
        0.8515 0.60432 0.05583
        0.23623 0.51395 0.69417
        0.49177 0.04136 0.4901
        0.08018 0.38669 0.18406
        0.78164 0.02659 0.94073
        0.58541 0.43914 0.14958
        0.42696 0.76378 0.43083
        0.73623 0.48605 0.05583
        0.80692 0.35523 0.99677
        0.76378 0.01395 0.55583
        0.35523 0.80692 0.00323
        0.9396 0.28164 0.56573
        0.9074 0.62573 0.80906
        0.51395 0.23623 0.30583
        0.3515 0.39568 0.69417
        0.47782 0.63159 0.12177
        0.43914 0.58541 0.85042
        0.68823 0.83545 0.97219
        0.78265 0.14377 0.09698
        0.47782 0.86841 0.12823
        0.70937 0.74055 0.9375
        0.92778 0.50241 0.35042
        0.18823 0.16455 0.77781
        0.41459 0.56086 0.14958
        0.97342 0.21836 0.05927
        0.20278 0.05755 0.97198
        0.1485 0.39568 0.05583
        0.00241 0.07223 0.10042
        0.8636 0.3636 0.25
        0.64791 0.09808 0.00667
        0.57223 0.50241 0.39958
        0.85209 0.40192 0.50667
        0.78164 0.4396 0.93427
        0.76378 0.42696 0.56917
        0.43433 0.64055 0.59719
        0.4396 0.71836 0.18427
        0.58018 0.61332 0.56594
        0.50823 0.54136 0.7599
        0.13951 0.29618 0.87156
        0.85791 0.49187 0.72219
        0.52659 0.21836 0.69073
        0.18823 0.33545 0.47219
        0.92881 0.77609 0.86833
        0.29063 0.24055 0.3125
        0.88669 0.58018 0.31594
        0.44246 0.29723 0.52802
        0.97782 0.36841 0.62823
        0.05755 0.29723 0.22198
        0.21736 0.64377 0.15302
        0.75946 0.70937 0.6875
        0.42778 0.4976 0.39958
        0.52218 0.13159 0.12823
        0.35945 0.93433 0.84719
        0.08541 0.56086 0.60042
        0.89841 0.55741 0.93083
        0.21836 0.52659 0.30927
        0.9074 0.87427 0.44094
        0.42881 0.27609 0.36833
        0.8636 0.1364 0.0
        0.88972 0.75 0.125
        0.64209 0.49187 0.02781
        0.45864 0.49177 0.2401
        0.66455 0.68823 0.72219
        0.44259 0.39841 0.18083
        0.99177 0.54136 0.9901
        0.58018 0.88669 0.68406
        0.22472 0.77327 0.3125
        0.61332 0.58018 0.43406
        0.98605 0.26378 0.80583
        0.56568 0.14055 0.65281
        0.52218 0.36841 0.12177
        0.87427 0.5926 0.69094
        0.33545 0.18823 0.52781
        0.62573 0.9074 0.19094
        0.43433 0.85945 0.65281
        0.55755 0.70278 0.52802
        0.25946 0.20937 0.1875
        0.73623 0.42696 0.18083
        0.09808 0.85209 0.25667
        0.29618 0.13951 0.12844
        0.14377 0.71736 0.34698
        0.79723 0.94246 0.97198
        0.14055 0.56568 0.34719
        0.51395 0.26378 0.94417
        0.54136 0.99177 0.0099
        0.94246 0.70278 0.22198
        0.93914 0.08541 0.35042
        0.44259 0.10159 0.06917
        0.99187 0.35791 0.77781
        0.79723 0.55755 0.27802
        0.3636 0.6364 0.5
        0.14209 0.99187 0.52781
        0.77609 0.57119 0.11833
        0.25946 0.29063 0.0625
        0.35624 0.71736 0.40302
        0.72472 0.27327 0.8125
        0.74055 0.79063 0.1875
        0.02659 0.78164 0.05927
        0.98605 0.23623 0.44417
        0.20937 0.25946 0.8125
        0.83545 0.68823 0.02781
        0.50241 0.57223 0.60042
        0.05755 0.20278 0.02802
        0.88669 0.91982 0.93406
        0.07223 0.00241 0.89958
        0.64377 0.21736 0.84698
        0.64055 0.06568 0.84719
        0.39568 0.3515 0.30583
        0.24055 0.20937 0.5625
        0.36841 0.97782 0.37177
        0.31177 0.16455 0.97219
        0.00823 0.45864 0.9901
        0.38669 0.08018 0.81594
        0.27327 0.77528 0.0625
        0.85624 0.28265 0.34698
        0.77528 0.22673 0.3125
        0.54136 0.50823 0.2401
        0.81177 0.83545 0.77781
        0.99177 0.95864 0.2599
        0.78265 0.35624 0.15302
        0.23623 0.98605 0.55583
        0.07304 0.76378 0.31917
        0.61028 0.75 0.625
        0.1364 0.6364 0.25
        0.38669 0.41982 0.43406
        0.85945 0.43433 0.34719
        0.07119 0.22391 0.86833
        0.81177 0.66455 0.47219
        0.4074 0.37427 0.94094
        0.95864 0.99177 0.7401
        0.87427 0.9074 0.55906
        0.50823 0.95864 0.4901
        0.20937 0.24055 0.4375
        0.40192 0.85209 0.49333
        0.91459 0.06086 0.64958
        0.57119 0.77609 0.88167
        0.31177 0.33545 0.27781
        0.91982 0.88669 0.06594
        0.91982 0.61332 0.18406
        0.48605 0.76378 0.30583
        0.11332 0.41982 0.31594
        0.47342 0.71836 0.55927
        0.63159 0.47782 0.87823
        0.93433 0.35945 0.15281
        0.45864 0.00823 0.0099
        0.12573 0.4074 0.69094
        0.57304 0.26378 0.81917
        0.0 0.5 0.25
        0.6364 0.1364 0.75
        0.72673 0.27528 0.1875
        0.28265 0.85624 0.65302
        0.79618 0.86049 0.62156
        0.19308 0.85523 0.25323
        0.05741 0.10159 0.68083
        0.35791 0.99187 0.22219
        0.04136 0.49177 0.5099
        0.50813 0.14209 0.27781
        0.05741 0.39841 0.56917
        0.57119 0.72391 0.36833
        0.1364 0.8636 0.0
        0.14477 0.80692 0.74677
        0.0 0.0 0.0
        0.89841 0.94259 0.31917
        0.78164 0.0604 0.31573
        0.72391 0.92881 0.61833
        0.75946 0.79063 0.5625
        0.90192 0.35209 0.99333
        0.69308 0.14477 0.49677
        0.64477 0.19308 0.00323
        0.33545 0.31177 0.72219
        0.5926 0.87427 0.30906
        0.92881 0.72391 0.38167
        0.24055 0.29063 0.6875
        0.28164 0.52659 0.44073
        0.64377 0.28265 0.40302
        0.13951 0.20382 0.37844
        0.26378 0.51395 0.05583
        0.30692 0.85523 0.49677
        0.06086 0.91459 0.35042
        0.69308 0.35523 0.75323
        0.79063 0.75946 0.4375
        """

        self.coord = "relative"

        self.cages = """
        12 0.03291 0.16545 0.26375
        14 0.0 0.0 0.5
        15 1.1529 0.6529 0.75
        14 0.87655 0.47671 1.11208
        16 0.75 0.91252 1.375
        12 0.66545 0.53291 0.23625
        12 0.33455 1.03291 1.01375
        14 0.62345 0.47671 0.63792
        16 0.08748 1.25 1.625
        12 0.53291 0.83455 0.48625
        14 0.12345 -0.02329 0.13792
        12 -0.75 0.64887 0.375
        14 0.62345 1.02329 0.61208
        14 -0.02329 1.12345 0.86208
        16 -0.08748 0.75 0.625
        12 1.16545 1.03291 0.73625
        14 0.37655 0.52329 0.63792
        12 0.96709 0.66545 0.98625
        12 -0.16545 0.53291 0.51375
        14 0.5 0.0 0.25
        14 -0.52329 0.62345 0.36208
        12 -0.18765 0.68765 0.25
        14 0.12345 0.52329 0.11208
        14 0.0 0.5 0.75
        12 0.53291 0.66545 0.76375
        12 0.46709 0.16545 0.48625
        16 0.75 0.58748 0.875
        12 0.14887 0.75 0.125
        14 0.5 0.5 1.0
        12 0.18765 0.31235 0.25
        14 0.02329 0.87655 0.86208
        15 0.3471 0.8471 1.25
        14 1.52329 1.12345 0.88792
        15 -0.3471 0.3471 1.0
        12 -0.16545 0.96709 0.73625
        12 0.03291 0.33455 -0.01375
        12 0.16545 1.46709 1.51375
        15 0.1529 0.8471 1.5
        12 0.35113 -0.25 0.625
        15 -0.1529 0.3471 0.75
        12 -0.25 1.14887 0.875
        12 0.33455 1.46709 1.23625
        14 1.02329 0.62345 0.38792
        16 0.41252 1.25 1.125
        12 0.64887 1.25 0.625
        14 -0.02329 1.37655 1.38792
        12 -0.33455 0.96709 1.01375
        16 -0.41252 0.75 1.125
        12 0.68765 0.68765 0.5
        15 0.6529 0.1529 0.25
        14 0.37655 -0.02329 0.61208
        16 0.25 0.41252 0.875
        12 0.81235 0.81235 1.0
        12 -0.31235 0.81235 0.75
        12 0.96709 0.83455 1.26375
        12 0.31235 0.31235 0.5
        15 -0.1529 1.1529 0.5
        12 0.18765 0.18765 0.0
        12 0.75 1.35113 1.375
        12 0.46709 0.33455 0.76375
        16 0.25 0.08748 0.375
        12 0.85113 0.25 1.125
        14 0.87655 1.02329 1.13792
        12 0.31235 0.18765 0.75
        14 0.52329 1.37655 1.36208
        15 0.3471 1.6529 1.0
        14 -0.52329 0.87655 0.88792
        12 0.25 0.85113 0.875
        """

        self.bondlen = 3

        self.cell = """
        32.72400861793071 32.72400861793071 19.43647231252596
        """

        self.density = 0.5572024775190785

        self.cell = cellvectors(a=32.72400861793071,
                                b=32.72400861793071,
                                c=19.43647231252596)
