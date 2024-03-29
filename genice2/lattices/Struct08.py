# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (40,20,8,12,)
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
        280 352
        417 212
        306 438
        110 27
        137 248
        52 289
        423 238
        170 433
        399 417
        82 168
        410 442
        208 199
        154 172
        46 294
        114 362
        69 351
        104 143
        180 12
        15 48
        13 323
        209 76
        242 162
        435 17
        151 251
        195 376
        44 13
        406 354
        8 379
        29 141
        242 199
        324 135
        359 116
        230 321
        107 271
        425 408
        279 89
        453 438
        339 160
        417 134
        246 356
        145 255
        97 151
        108 392
        331 48
        207 321
        198 141
        159 334
        129 344
        136 222
        355 112
        122 357
        157 317
        347 148
        277 224
        246 420
        388 305
        182 297
        234 113
        54 270
        0 134
        425 241
        96 157
        240 402
        390 144
        282 201
        94 71
        90 247
        450 408
        227 26
        282 361
        259 311
        292 323
        299 308
        369 155
        119 454
        276 201
        226 24
        165 358
        275 432
        215 102
        6 345
        149 321
        42 267
        143 439
        77 371
        203 400
        27 92
        163 161
        446 411
        57 351
        137 309
        354 60
        131 123
        385 73
        264 372
        88 168
        286 199
        176 28
        453 274
        322 256
        431 190
        103 423
        429 288
        9 434
        150 325
        40 146
        292 300
        81 33
        24 244
        418 344
        68 258
        343 55
        99 121
        151 268
        237 149
        357 268
        204 363
        223 221
        118 55
        35 8
        20 332
        62 454
        214 221
        114 11
        43 314
        85 338
        249 176
        367 342
        219 124
        150 168
        368 426
        325 60
        300 449
        324 448
        192 2
        80 134
        261 123
        166 443
        173 407
        453 217
        202 391
        264 90
        213 352
        406 380
        341 101
        47 455
        110 101
        104 236
        340 254
        241 125
        397 139
        356 14
        403 394
        233 220
        399 3
        234 18
        203 188
        398 379
        387 91
        233 180
        233 225
        444 436
        336 322
        158 131
        45 73
        265 130
        117 328
        322 23
        410 263
        309 413
        403 15
        260 227
        46 350
        160 288
        9 84
        163 361
        108 285
        375 448
        44 361
        428 405
        103 390
        231 149
        191 212
        130 418
        132 46
        133 248
        291 144
        392 39
        132 139
        282 211
        41 59
        435 247
        429 389
        318 394
        167 138
        374 452
        293 269
        111 36
        200 313
        87 243
        132 332
        155 308
        413 224
        86 287
        132 235
        443 164
        160 227
        440 25
        375 79
        343 245
        274 338
        420 232
        206 53
        417 348
        302 298
        152 400
        398 189
        45 4
        286 0
        208 267
        10 288
        17 445
        147 21
        367 262
        281 327
        158 66
        319 373
        291 437
        119 38
        3 341
        447 441
        305 73
        402 68
        201 415
        205 61
        418 25
        322 267
        145 198
        356 18
        340 370
        47 176
        83 54
        45 5
        397 350
        265 386
        31 307
        266 20
        207 142
        126 363
        0 430
        62 143
        99 47
        82 9
        31 142
        163 409
        3 293
        158 197
        140 63
        36 289
        379 14
        310 80
        215 414
        264 299
        325 78
        157 216
        44 315
        169 49
        452 304
        59 394
        103 429
        81 12
        401 151
        148 323
        439 381
        138 67
        29 318
        165 16
        281 295
        17 371
        314 65
        56 324
        398 390
        336 452
        307 21
        409 29
        62 329
        291 189
        276 421
        7 182
        369 365
        209 220
        283 37
        143 313
        377 247
        229 257
        275 454
        98 252
        74 329
        419 307
        435 447
        111 58
        50 363
        94 75
        202 138
        425 352
        260 225
        136 334
        121 428
        340 197
        195 248
        238 18
        272 320
        392 326
        437 455
        101 333
        192 414
        99 439
        40 405
        6 224
        362 22
        309 371
        415 221
        253 192
        378 447
        75 61
        218 159
        303 397
        41 317
        393 420
        71 365
        49 210
        122 290
        346 353
        201 126
        1 319
        50 4
        374 348
        83 345
        283 212
        43 365
        173 332
        191 220
        264 302
        194 307
        102 135
        399 427
        285 2
        396 128
        205 382
        295 84
        296 85
        396 145
        88 70
        171 353
        286 416
        112 13
        431 238
        409 158
        41 370
        125 129
        39 92
        87 380
        337 262
        301 73
        284 68
        337 54
        219 363
        50 100
        166 427
        1 317
        335 323
        16 441
        100 425
        284 183
        343 222
        365 120
        93 368
        53 115
        101 217
        153 281
        130 380
        251 426
        227 326
        205 279
        331 45
        1 397
        352 344
        223 386
        253 319
        166 180
        189 80
        64 364
        412 414
        419 424
        107 383
        153 308
        121 239
        38 432
        273 169
        6 97
        286 182
        451 453
        95 358
        280 25
        157 394
        366 129
        130 183
        242 256
        26 81
        67 178
        171 179
        5 444
        276 335
        159 376
        349 161
        369 279
        427 113
        50 48
        277 239
        43 364
        55 355
        194 385
        152 326
        108 27
        69 267
        186 321
        213 450
        288 330
        410 224
        297 134
        91 9
        243 159
        387 106
        370 289
        330 170
        203 350
        196 181
        65 206
        372 172
        368 328
        252 241
        34 13
        318 66
        71 116
        137 441
        348 220
        327 154
        93 80
        110 293
        84 154
        265 354
        272 401
        229 376
        181 179
        88 222
        8 28
        98 169
        393 304
        401 337
        31 72
        218 422
        273 366
        10 423
        233 274
        47 144
        231 48
        395 28
        281 360
        7 336
        76 452
        131 177
        11 216
        7 368
        253 41
        166 348
        165 383
        446 2
        338 235
        160 333
        139 39
        94 313
        211 29
        115 351
        377 372
        44 181
        153 302
        96 20
        393 138
        347 168
        97 226
        86 341
        207 5
        272 53
        67 416
        209 287
        229 133
        147 261
        277 382
        381 200
        337 271
        455 268
        81 330
        253 32
        2 278
        205 71
        442 378
        339 341
        271 69
        399 287
        207 156
        312 8
        56 444
        153 422
        440 30
        427 296
        377 263
        381 38
        291 244
        214 169
        404 206
        406 30
        89 17
        369 172
        316 426
        156 102
        401 345
        180 85
        117 312
        367 53
        245 315
        15 126
        447 263
        19 74
        300 112
        96 141
        415 353
        395 238
        193 125
        147 177
        450 21
        378 445
        37 269
        255 161
        195 78
        64 320
        251 24
        186 311
        156 72
        332 319
        162 342
        240 243
        93 437
        292 250
        188 296
        19 313
        416 190
        285 373
        118 125
        72 128
        403 436
        120 372
        408 123
        249 432
        358 314
        328 297
        357 262
        52 414
        423 420
        226 38
        214 124
        174 148
        119 146
        7 242
        204 252
        95 61
        171 354
        15 421
        371 133
        87 440
        272 162
        386 187
        245 196
        152 225
        213 300
        92 51
        211 210
        259 173
        197 185
        375 412
        311 216
        121 268
        339 312
        240 55
        109 67
        276 174
        360 133
        72 185
        108 333
        174 221
        404 262
        198 59
        111 22
        75 57
        301 79
        129 380
        186 301
        66 254
        424 388
        451 411
        422 441
        39 105
        109 104
        37 217
        355 315
        161 21
        396 52
        246 443
        92 36
        273 406
        345 206
        74 351
        239 381
        113 225
        122 226
        416 164
        297 287
        395 391
        119 428
        230 444
        419 255
        177 315
        349 66
        105 235
        219 366
        1 11
        79 22
        178 199
        58 192
        20 32
        384 407
        31 261
        266 139
        208 342
        152 296
        432 244
        228 357
        26 184
        203 451
        94 364
        343 449
        183 360
        429 312
        188 294
        74 42
        200 116
        335 361
        110 451
        385 185
        377 77
        209 113
        306 278
        237 126
        100 98
        405 382
        107 364
        123 185
        64 115
        37 170
        10 12
        292 245
        58 370
        409 145
        204 149
        100 124
        117 310
        215 388
        384 79
        104 23
        90 89
        223 347
        107 57
        347 34
        96 259
        147 193
        439 249
        278 412
        246 12
        43 140
        407 448
        305 135
        254 255
        230 127
        410 359
        284 229
        260 18
        109 391
        336 430
        58 388
        290 442
        83 63
        34 179
        448 306
        299 279
        263 309
        202 232
        389 175
        76 190
        283 3
        269 438
        33 105
        106 34
        339 184
        431 304
        165 308
        421 210
        265 295
        421 318
        400 303
        316 176
        97 239
        280 449
        236 256
        11 56
        16 270
        216 436
        385 424
        311 384
        283 175
        40 61
        434 172
        373 235
        390 232
        78 154
        194 408
        359 40
        103 28
        298 60
        383 140
        65 63
        69 367
        19 236
        14 232
        359 120
        90 91
        247 334
        248 327
        184 175
        64 19
        109 431
        278 51
        77 434
        211 237
        218 435
        258 257
        254 59
        310 212
        325 30
        304 23
        51 269
        273 353
        374 164
        189 35
        280 118
        275 202
        306 411
        294 338
        6 290
        383 54
        26 356
        362 324
        445 413
        241 261
        214 187
        375 102
        398 389
        316 244
        83 413
        95 299
        175 330
        231 210
        204 142
        222 258
        183 243
        303 392
        76 182
        346 282
        194 4
        131 250
        188 293
        63 137
        422 195
        442 405
        75 329
        167 42
        135 5
        440 376
        106 82
        187 30
        234 443
        122 428
        271 320
        275 144
        266 289
        302 78
        99 454
        305 22
        404 270
        200 146
        326 33
        266 317
        33 85
        184 379
        402 112
        396 197
        196 148
        0 374
        362 411
        86 191
        124 344
        142 4
        349 335
        450 250
        70 60
        89 434
        331 436
        36 412
        346 49
        106 402
        251 162
        438 294
        407 46
        395 14
        433 105
        236 178
        52 32
        146 329
        150 84
        257 77
        215 128
        228 437
        93 430
        27 446
        42 178
        257 334
        173 56
        218 360
        127 141
        16 65
        95 120
        314 115
        98 231
        228 24
        62 167
        127 237
        35 328
        174 49
        426 455
        349 250
        400 274
        252 366
        198 32
        234 190
        140 155
        285 433
        196 88
        228 342
        136 387
        277 116
        404 290
        358 57
        208 430
        449 68
        191 217
        114 384
        170 333
        223 171
        340 424
        163 177
        433 51
        156 301
        445 382
        260 10
        387 70
        219 415
        259 230
        111 446
        70 179
        378 270
        391 249
        303 373
        82 258
        193 355
        418 187
        213 193
        403 127
        155 327
        331 186
        136 240
        298 295
        114 350
        256 320
        284 25
        167 23
        298 91
        316 35
        117 86
        118 87
        389 310
        393 164
        150 386
        419 128
        346 181
        """

        self.waters = """
        0.30758 0.06551 0.11251
        0.81035 0.32455 0.8286
        0.0 0.29662 0.47427
        0.5 0.18246 0.25733
        0.875 0.47084 0.19269
        0.0 0.40723 0.13518
        0.625 0.88098 0.38784
        0.81742 0.03098 0.14926
        0.69158 0.08983 0.46528
        0.875 0.71215 0.72358
        0.0 0.14663 0.70953
        0.80843 0.34212 0.95286
        0.18343 0.15808 0.76208
        0.0 0.59277 0.63518
        0.5 0.09277 0.63518
        0.80758 0.47231 0.91257
        0.30843 0.81886 0.16015
        0.0 0.77685 0.48408
        0.69158 0.12273 0.77692
        0.625 0.93411 0.87511
        0.31466 0.35473 0.75205
        0.0 0.50507 0.4329
        0.68343 0.34192 0.26208
        0.0 0.00549 0.88518
        0.5 0.97519 0.3718
        0.68535 0.64527 0.25205
        0.5 0.15913 0.62884
        0.68535 0.26342 0.42177
        0.81658 0.06553 0.5584
        0.375 0.46061 0.75393
        0.5 0.65147 0.07264
        0.19242 0.4755 0.30246
        0.18966 0.36899 0.66596
        0.375 0.21215 0.72358
        0.0 0.6239 0.70549
        0.625 0.06589 0.37511
        0.5 0.31824 0.46101
        0.18535 0.20691 0.33748
        0.5 0.95055 0.56158
        0.5 0.26776 0.64121
        0.31035 0.86899 0.66596
        0.81035 0.36899 0.66596
        0.30758 0.97231 0.91257
        0.69158 0.84212 0.95286
        0.19158 0.56812 0.65786
        0.80843 0.43188 0.15786
        0.375 0.28865 0.99304
        0.0 0.0081 0.51018
        0.68258 0.47302 0.02541
        0.5 0.55128 0.88493
        0.80758 0.49605 0.11357
        0.31466 0.26342 0.42177
        0.31466 0.36981 0.55783
        0.5 0.90723 0.13518
        0.0 0.85338 0.20953
        0.31658 0.62006 0.43005
        0.0 0.35268 0.00585
        0.18343 0.87995 0.93005
        0.81035 0.35007 0.46744
        0.875 0.4105 0.70005
        0.31466 0.67565 0.88782
        0.18535 0.85473 0.75205
        0.125 0.96061 0.75393
        0.69158 0.81886 0.16015
        0.69158 0.91018 0.96528
        0.5 0.84087 0.12884
        0.68258 0.46902 0.64926
        0.5 0.02481 0.8718
        0.80843 0.65788 0.45286
        0.18343 0.93447 0.0584
        0.31658 0.65808 0.76208
        0.81466 0.85473 0.75205
        0.31742 0.43772 0.33906
        0.68343 0.41338 0.24704
        0.375 0.93411 0.87511
        0.125 0.88823 0.82299
        0.81742 0.08828 0.0241
        0.68966 0.74534 0.54677
        0.5 0.72315 0.98408
        0.5 0.35338 0.20953
        0.25 0.07728 0.28581
        0.30843 0.18115 0.66015
        0.80843 0.68115 0.66015
        0.81658 0.84192 0.26208
        0.81035 0.70783 0.87983
        0.31035 0.20783 0.87983
        0.81466 0.14527 0.25205
        0.31466 0.64527 0.25205
        0.5 0.64663 0.70953
        0.0 0.76776 0.64121
        0.18535 0.74524 0.67177
        0.125 0.71215 0.72358
        0.5 0.27685 0.48408
        0.125 0.03939 0.25393
        0.875 0.88823 0.82299
        0.31035 0.82455 0.8286
        0.375 0.38823 0.82299
        0.69242 0.91678 0.4273
        0.5 0.53058 0.08144
        0.0 0.96942 0.58144
        0.68258 0.53098 0.14926
        0.81466 0.20691 0.33748
        0.19158 0.37727 0.27692
        0.0 0.0852 0.57616
        0.81742 0.99088 0.828
        0.31466 0.24524 0.67177
        0.0 0.65913 0.62884
        0.0 0.87959 0.00049
        0.81035 0.24534 0.54677
        0.69242 0.0245 0.80246
        0.68966 0.23645 0.29677
        0.68535 0.32435 0.38782
        0.0 0.60575 0.5186
        0.68966 0.14993 0.96744
        0.68343 0.31976 0.05461
        0.5 0.89425 0.0186
        0.68966 0.86899 0.66596
        0.875 0.11177 0.32299
        0.375 0.61177 0.32299
        0.30758 0.93449 0.61251
        0.5 0.81754 0.75733
        0.0 0.93235 0.50519
        0.30758 0.91678 0.4273
        0.5 0.49451 0.38518
        0.80758 0.56551 0.11251
        0.25 0.57728 0.28581
        0.0 0.49493 0.9329
        0.125 0.43411 0.87511
        0.19242 0.41678 0.4273
        0.125 0.5895 0.20005
        0.0 0.65686 0.14876
        0.5 0.5081 0.51018
        0.31466 0.2931 0.83748
        0.81035 0.73645 0.29677
        0.375 0.0895 0.20005
        0.0 0.3761 0.20549
        0.31658 0.68025 0.55461
        0.625 0.78785 0.22358
        0.30758 0.0245 0.80246
        0.5 0.29742 0.74373
        0.81658 0.81976 0.05461
        0.25 0.42272 0.78581
        0.125 0.47084 0.19269
        0.875 0.96061 0.75393
        0.18258 0.02698 0.52541
        0.19242 0.43449 0.61251
        0.375 0.9105 0.70005
        0.19242 0.52769 0.41257
        0.68343 0.58662 0.74704
        0.31742 0.47302 0.02541
        0.68535 0.67565 0.88782
        0.81742 0.93772 0.33906
        0.68966 0.20783 0.87983
        0.18966 0.75467 0.04677
        0.68535 0.73659 0.92177
        0.875 0.78865 0.99304
        0.31658 0.41338 0.24704
        0.625 0.38823 0.82299
        0.5 0.46942 0.58144
        0.31466 0.70691 0.33748
        0.81658 0.18025 0.55461
        0.0 0.49005 0.5579
        0.625 0.97084 0.19269
        0.19242 0.50395 0.61357
        0.30758 0.08323 0.9273
        0.18343 0.81976 0.05461
        0.31035 0.14993 0.96744
        0.18258 0.99088 0.828
        0.68343 0.65808 0.76208
        0.5 0.56765 0.00519
        0.125 0.21135 0.49304
        0.125 0.61902 0.88784
        0.68966 0.76356 0.79677
        0.19158 0.34212 0.95286
        0.68258 0.56228 0.83906
        0.30843 0.15788 0.45286
        0.81742 0.02698 0.52541
        0.31742 0.52698 0.52541
        0.5 0.99493 0.9329
        0.19158 0.62273 0.77692
        0.18535 0.17565 0.88782
        0.31658 0.58662 0.74704
        0.69242 0.06551 0.11251
        0.0 0.68246 0.25733
        0.5 0.14732 0.50585
        0.5 0.44872 0.38493
        0.5 0.4148 0.07616
        0.68535 0.63019 0.05783
        0.5 0.23202 0.07802
        0.375 0.06589 0.37511
        0.69242 0.08323 0.9273
        0.0 0.16319 0.19823
        0.0 0.33133 0.52175
        0.125 0.56589 0.37511
        0.80758 0.4755 0.30246
        0.5 0.73224 0.14121
        0.5 0.59996 0.7009
        0.5 0.43235 0.50519
        0.125 0.4105 0.70005
        0.5 0.00995 0.0579
        0.625 0.9105 0.70005
        0.0 0.52481 0.8718
        0.375 0.02917 0.69269
        0.68966 0.25467 0.04677
        0.19242 0.49605 0.11357
        0.0 0.83681 0.69823
        0.5 0.8761 0.20549
        0.19158 0.43188 0.15786
        0.30758 0.99605 0.11357
        0.81466 0.13019 0.05783
        0.5 0.50549 0.88518
        0.31742 0.49088 0.828
        0.18535 0.14527 0.25205
        0.875 0.56589 0.37511
        0.68258 0.58828 0.0241
        0.125 0.38098 0.38784
        0.68343 0.37995 0.93005
        0.0 0.20258 0.24373
        0.18966 0.73645 0.29677
        0.0 0.54946 0.06158
        0.0 0.15147 0.07264
        0.80758 0.58323 0.9273
        0.5 0.65835 0.58512
        0.875 0.61902 0.88784
        0.68966 0.85007 0.46744
        0.81466 0.17565 0.88782
        0.5 0.93995 0.44789
        0.69158 0.18115 0.66015
        0.30758 0.9755 0.30246
        0.68535 0.70691 0.33748
        0.19158 0.41018 0.96528
        0.5 0.4919 0.01018
        0.30843 0.06812 0.65786
        0.0 0.18177 0.96101
        0.625 0.11902 0.88784
        0.18966 0.26356 0.79677
        0.69242 0.97231 0.91257
        0.19242 0.47231 0.91257
        0.81658 0.08662 0.74704
        0.81742 0.91172 0.5241
        0.19158 0.65788 0.45286
        0.375 0.53939 0.25393
        0.69242 0.99605 0.11357
        0.18966 0.67545 0.3286
        0.5 0.00507 0.4329
        0.5 0.5852 0.57616
        0.30843 0.12273 0.77692
        0.31035 0.74534 0.54677
        0.68535 0.75477 0.17177
        0.69242 0.00395 0.61357
        0.68258 0.52698 0.52541
        0.69242 0.9755 0.30246
        0.31742 0.53098 0.14926
        0.0 0.34314 0.64876
        0.80758 0.43449 0.61251
        0.0 0.45055 0.56158
        0.81742 0.97302 0.02541
        0.625 0.71135 0.49304
        0.68343 0.68025 0.55461
        0.31658 0.37995 0.93005
        0.81658 0.15808 0.76208
        0.31742 0.50913 0.328
        0.18343 0.91338 0.24704
        0.5 0.79662 0.47427
        0.31035 0.76356 0.79677
        0.0 0.66867 0.02175
        0.5 0.33681 0.69823
        0.18258 0.97302 0.02541
        0.0 0.94872 0.38493
        0.31035 0.23645 0.29677
        0.18343 0.84192 0.26208
        0.0 0.9148 0.07616
        0.69158 0.93188 0.15786
        0.31742 0.58828 0.0241
        0.0 0.22315 0.98408
        0.30758 0.00395 0.61357
        0.80758 0.5245 0.80246
        0.81466 0.86981 0.55783
        0.18966 0.29218 0.37983
        0.0 0.79742 0.74373
        0.625 0.61177 0.32299
        0.0 0.73202 0.07802
        0.19242 0.5245 0.80246
        0.31035 0.17545 0.3286
        0.81035 0.67545 0.3286
        0.0 0.26798 0.57802
        0.5 0.04946 0.06158
        0.68966 0.13101 0.16596
        0.0 0.15835 0.58512
        0.5 0.34853 0.57264
        0.375 0.88098 0.38784
        0.30758 0.02769 0.41257
        0.68343 0.56553 0.5584
        0.5 0.22291 0.23408
        0.31035 0.25467 0.04677
        0.0 0.70338 0.97427
        0.5 0.20338 0.97427
        0.625 0.0895 0.20005
        0.18966 0.70783 0.87983
        0.18535 0.7931 0.83748
        0.80843 0.58983 0.46528
        0.5 0.40004 0.2009
        0.31466 0.73659 0.92177
        0.81035 0.26356 0.79677
        0.0 0.05128 0.88493
        0.80843 0.37727 0.27692
        0.125 0.28785 0.22358
        0.0 0.47519 0.3718
        0.125 0.78865 0.99304
        0.68966 0.79218 0.37983
        0.125 0.11177 0.32299
        0.5 0.37959 0.00049
        0.81658 0.12006 0.43005
        0.75 0.92272 0.78581
        0.5 0.85268 0.00585
        0.31658 0.56553 0.5584
        0.69242 0.02769 0.41257
        0.68535 0.35473 0.75205
        0.625 0.46061 0.75393
        0.0 0.31754 0.75733
        0.81658 0.93447 0.0584
        0.31658 0.43447 0.0584
        0.0 0.9919 0.01018
        0.80843 0.56812 0.65786
        0.0 0.34087 0.12884
        0.5 0.68177 0.96101
        0.625 0.21215 0.72358
        0.81035 0.75467 0.04677
        0.75 0.07728 0.28581
        0.25 0.92272 0.78581
        0.18343 0.18025 0.55461
        0.68343 0.43447 0.0584
        0.18966 0.32455 0.8286
        0.875 0.21135 0.49304
        0.375 0.71135 0.49304
        0.875 0.52917 0.69269
        0.0 0.03058 0.08144
        0.0 0.90004 0.2009
        0.18535 0.23659 0.92177
        0.69158 0.15788 0.45286
        0.68258 0.41172 0.5241
        0.68966 0.17545 0.3286
        0.375 0.97084 0.19269
        0.5 0.62041 0.50049
        0.875 0.5895 0.20005
        0.69158 0.87727 0.27692
        0.31742 0.56228 0.83906
        0.80843 0.62273 0.77692
        0.18535 0.13019 0.05783
        0.80758 0.50395 0.61357
        0.625 0.28865 0.99304
        0.30843 0.91018 0.96528
        0.75 0.57728 0.28581
        0.19242 0.58323 0.9273
        0.18966 0.64993 0.96744
        0.19158 0.58983 0.46528
        0.5 0.1239 0.70549
        0.18258 0.93772 0.33906
        0.30843 0.84212 0.95286
        0.5 0.84314 0.64876
        0.0 0.72291 0.23408
        0.125 0.52917 0.69269
        0.80843 0.31886 0.16015
        0.0 0.50995 0.0579
        0.81658 0.87995 0.93005
        0.68966 0.82455 0.8286
        0.19242 0.56551 0.11251
        0.30843 0.93188 0.15786
        0.875 0.03939 0.25393
        0.81466 0.7931 0.83748
        0.68535 0.36981 0.55783
        0.81466 0.76342 0.42177
        0.5 0.77709 0.73408
        0.0 0.27709 0.73408
        0.18258 0.08828 0.0241
        0.31658 0.34192 0.26208
        0.5 0.70258 0.24373
        0.5 0.76798 0.57802
        0.18535 0.82435 0.38782
        0.5 0.10575 0.5186
        0.18966 0.63101 0.16596
        0.69242 0.93449 0.61251
        0.0 0.84853 0.57264
        0.0 0.84165 0.08512
        0.5 0.34165 0.08512
        0.68258 0.43772 0.33906
        0.81035 0.64993 0.96744
        0.19158 0.68115 0.66015
        0.875 0.38098 0.38784
        0.18343 0.12006 0.43005
        0.18343 0.06553 0.5584
        0.625 0.02917 0.69269
        0.68535 0.24524 0.67177
        0.18258 0.06228 0.83906
        0.75 0.42272 0.78581
        0.69158 0.06812 0.65786
        0.31742 0.41172 0.5241
        0.68535 0.2931 0.83748
        0.30843 0.08983 0.46528
        0.5 0.15686 0.14876
        0.81466 0.23659 0.92177
        0.81658 0.91338 0.24704
        0.0 0.64732 0.50585
        0.875 0.43411 0.87511
        0.30843 0.87727 0.27692
        0.18535 0.86981 0.55783
        0.31466 0.63019 0.05783
        0.31658 0.31976 0.05461
        0.68258 0.50913 0.328
        0.31742 0.46902 0.64926
        0.5 0.83133 0.52175
        0.875 0.28785 0.22358
        0.31466 0.32435 0.38782
        0.81466 0.82435 0.38782
        0.18966 0.35007 0.46744
        0.0 0.56006 0.94789
        0.5 0.06006 0.94789
        0.31035 0.13101 0.16596
        0.81035 0.63101 0.16596
        0.0 0.43995 0.44789
        0.18343 0.08662 0.74704
        0.68258 0.49088 0.828
        0.31466 0.75477 0.17177
        0.0 0.09996 0.7009
        0.80758 0.41678 0.4273
        0.625 0.53939 0.25393
        0.81742 0.00913 0.328
        0.5 0.16867 0.02175
        0.18258 0.91172 0.5241
        0.0 0.12041 0.50049
        0.18258 0.03098 0.14926
        0.81742 0.06228 0.83906
        0.5 0.99005 0.5579
        0.18966 0.24534 0.54677
        0.81466 0.74524 0.67177
        0.18535 0.76342 0.42177
        0.80843 0.41018 0.96528
        0.18258 0.00913 0.328
        0.18535 0.25477 0.17177
        0.81742 0.96902 0.64926
        0.5 0.66319 0.19823
        0.375 0.78785 0.22358
        0.31035 0.85007 0.46744
        0.375 0.11902 0.88784
        0.0 0.39425 0.0186
        0.0 0.81824 0.46101
        0.81035 0.29218 0.37983
        0.31035 0.79218 0.37983
        0.19158 0.31886 0.16015
        0.68343 0.62006 0.43005
        0.80758 0.52769 0.41257
        0.81466 0.25477 0.17177
        0.0 0.06765 0.00519
        0.0 0.23224 0.14121
        0.18258 0.96902 0.64926
        0.0 0.99451 0.38518
        """

        self.coord = "relative"

        self.cages = """
        14 0.5 0.12437 0.32504
        14 0.5 0.87563 0.82504
        14 0.0 0.62437 0.32504
        14 0.0 0.37563 0.82504
        14 0.0 0.06207 0.39351
        14 0.0 0.93793 0.89351
        14 0.5 0.56207 0.39351
        14 0.5 0.43793 0.89351
        12 0.0 0.45292 0.03577
        12 0.0 0.54708 0.53577
        12 0.5 0.95292 0.03577
        12 0.5 0.04708 0.53577
        16 0.0 0.02718 0.69751
        16 0.0 0.97282 0.19751
        16 0.5 0.52718 0.69751
        16 0.5 0.47282 0.19751
        12 0.0 0.2502 0.36055
        12 0.0 0.7498 0.86055
        12 0.5 0.7502 0.36055
        12 0.5 0.2498 0.86055
        14 0.0 0.09265 0.17855
        14 0.0 0.90735 0.67855
        14 0.5 0.59265 0.17855
        14 0.5 0.40735 0.67855
        15 0.5 0.3789 0.3925
        15 0.5 0.6211 0.8925
        15 0.0 0.8789 0.3925
        15 0.0 0.1211 0.8925
        12 0.0 0.42757 0.29136
        12 0.0 0.57243 0.79136
        12 0.5 0.92757 0.29136
        12 0.5 0.07243 0.79136
        12 0.5 0.11095 0.05049
        12 0.5 0.88905 0.55049
        12 0.0 0.61095 0.05049
        12 0.0 0.38905 0.55049
        16 0.5 0.21336 0.48476
        16 0.5 0.78664 0.98476
        16 0.0 0.71336 0.48476
        16 0.0 0.28664 0.98476
        12 0.26968 0.47158 0.47485
        12 0.73032 0.52842 0.97485
        12 0.26968 0.52842 0.97485
        12 0.76968 0.97158 0.47485
        12 0.73032 0.47158 0.47485
        12 0.23032 0.02842 0.97485
        12 0.76968 0.02842 0.97485
        12 0.23032 0.97158 0.47485
        15 0.5 0.03003 0.24613
        15 0.5 0.96997 0.74613
        15 0.0 0.53003 0.24613
        15 0.0 0.46997 0.74613
        12 0.0 0.168 0.39484
        12 0.0 0.832 0.89484
        12 0.5 0.668 0.39484
        12 0.5 0.332 0.89484
        16 0.5 0.28616 0.23202
        16 0.5 0.71384 0.73202
        16 0.0 0.78616 0.23202
        16 0.0 0.21384 0.73202
        12 0.2663 0.37422 0.1068
        12 0.7337 0.62578 0.6068
        12 0.2663 0.62578 0.6068
        12 0.7663 0.87422 0.1068
        12 0.7337 0.37422 0.1068
        12 0.2337 0.12578 0.6068
        12 0.7663 0.12578 0.6068
        12 0.2337 0.87422 0.1068
        14 0.24138 0.19606 0.10976
        14 0.75862 0.80394 0.60976
        14 0.24138 0.80394 0.60976
        14 0.74138 0.69606 0.10976
        14 0.75862 0.19606 0.10976
        14 0.25862 0.30394 0.60976
        14 0.74138 0.30394 0.60976
        14 0.25862 0.69606 0.10976
        12 0.0 0.3284 0.317
        12 0.0 0.6716 0.817
        12 0.5 0.8284 0.317
        12 0.5 0.1716 0.817
        """

        self.bondlen = 3

        self.cell = """
        14.682222842138556 80.57904303333389 24.251803115084247
        """

        self.density = 0.4750492601557878

        self.cell = cellvectors(a=14.682222842138556,
                                b=80.57904303333389,
                                c=24.251803115084247)
