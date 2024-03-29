# coding: utf-8
"""
Data source: Dutour Sikirić, Mathieu, Olaf Delgado-Friedrichs, and Michel Deza. “Space Fullerenes: a Computer Search for New Frank-Kasper Structures” Acta Crystallographica Section A Foundations of Crystallography 66.Pt 5 (2010): 602–615.

Cage composition:
 (12,14,15,16) = (40,16,16,8,)
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
        208 115
        328 222
        349 102
        337 211
        204 111
        84 288
        12 331
        171 314
        339 145
        70 174
        243 53
        350 91
        363 162
        335 190
        422 282
        423 283
        17 218
        142 440
        67 336
        366 425
        379 182
        403 284
        135 115
        100 93
        85 322
        205 454
        130 430
        170 399
        387 218
        383 14
        393 356
        454 336
        116 352
        311 410
        279 324
        288 97
        403 321
        264 79
        114 421
        276 294
        234 316
        220 191
        230 296
        376 166
        335 34
        23 175
        74 15
        432 441
        167 389
        173 262
        224 364
        154 140
        302 44
        373 192
        18 191
        154 69
        238 276
        33 38
        245 409
        78 450
        176 183
        129 75
        158 165
        62 63
        313 319
        97 438
        63 187
        269 359
        304 402
        44 48
        220 264
        135 18
        303 244
        178 375
        447 79
        395 80
        241 341
        26 213
        297 441
        249 252
        377 144
        419 183
        386 31
        269 335
        374 442
        157 396
        328 2
        189 89
        137 378
        327 24
        344 65
        171 232
        42 28
        315 306
        254 98
        203 392
        366 344
        266 79
        116 136
        95 246
        195 289
        120 228
        440 156
        29 401
        130 20
        259 164
        228 229
        203 292
        424 244
        408 179
        76 376
        163 429
        75 376
        117 361
        349 117
        80 382
        385 278
        159 67
        209 23
        82 384
        369 273
        199 175
        188 414
        141 320
        210 47
        285 227
        45 25
        43 277
        57 370
        86 279
        314 340
        254 95
        331 325
        241 400
        391 259
        194 58
        127 260
        257 176
        442 346
        228 71
        35 439
        359 41
        66 205
        190 27
        62 285
        254 421
        364 357
        354 238
        220 117
        429 36
        388 338
        138 193
        417 212
        87 362
        241 380
        220 354
        430 4
        433 217
        126 312
        200 103
        298 405
        19 417
        297 227
        12 277
        402 337
        127 131
        259 102
        11 424
        352 365
        372 193
        435 101
        444 305
        198 150
        437 252
        295 32
        27 261
        437 279
        418 109
        46 214
        399 420
        398 96
        270 441
        107 199
        161 252
        303 133
        412 324
        81 256
        158 110
        122 42
        59 275
        436 5
        13 31
        361 266
        170 143
        124 232
        422 8
        354 87
        13 282
        283 15
        171 9
        24 438
        448 405
        327 61
        58 389
        201 336
        161 289
        188 396
        443 201
        185 265
        139 451
        452 246
        59 50
        375 88
        347 353
        349 355
        294 206
        155 73
        13 118
        438 88
        212 101
        173 426
        366 330
        302 451
        268 189
        136 98
        176 86
        447 73
        370 407
        225 21
        219 310
        326 253
        209 166
        380 434
        210 166
        108 397
        105 55
        81 261
        433 405
        348 94
        151 449
        368 394
        225 155
        427 298
        192 78
        418 233
        309 426
        324 434
        99 89
        425 414
        448 454
        116 119
        109 47
        53 97
        143 21
        303 167
        367 326
        411 433
        178 378
        430 338
        108 256
        318 450
        54 336
        326 433
        122 405
        151 22
        115 58
        91 38
        52 110
        92 39
        219 308
        381 371
        382 399
        16 408
        37 292
        281 123
        386 239
        12 443
        82 340
        55 299
        60 300
        301 71
        181 216
        86 318
        335 262
        200 84
        160 146
        204 51
        321 272
        423 449
        70 316
        453 334
        122 237
        20 271
        213 234
        60 394
        284 334
        267 14
        74 148
        66 245
        11 371
        350 455
        40 17
        360 409
        258 156
        243 110
        423 340
        317 127
        181 182
        221 407
        331 394
        45 437
        221 400
        287 10
        349 55
        429 89
        138 147
        216 103
        107 304
        34 291
        25 145
        55 182
        198 92
        135 305
        22 445
        140 338
        451 374
        253 152
        342 190
        219 260
        449 211
        343 110
        125 353
        111 341
        271 47
        189 383
        68 326
        224 247
        453 410
        327 69
        307 92
        37 96
        199 390
        43 194
        416 10
        221 418
        410 365
        346 308
        128 87
        83 197
        403 265
        198 432
        263 172
        237 296
        171 2
        51 109
        46 234
        129 151
        269 297
        289 434
        257 450
        124 385
        178 93
        302 290
        250 438
        404 406
        401 399
        339 383
        195 233
        1 263
        427 253
        106 179
        20 140
        11 158
        210 57
        106 31
        95 192
        305 447
        401 127
        379 288
        369 119
        244 165
        0 251
        312 394
        294 115
        206 389
        180 307
        311 185
        191 371
        81 270
        242 39
        310 350
        355 329
        219 177
        369 334
        155 5
        232 445
        447 231
        80 8
        128 424
        388 212
        0 260
        422 169
        284 416
        28 346
        360 145
        268 90
        256 439
        54 49
        247 186
        124 116
        128 444
        417 446
        239 382
        262 180
        56 343
        0 267
        391 105
        386 332
        262 215
        257 313
        452 4
        314 356
        30 154
        351 319
        372 44
        435 70
        77 442
        338 446
        114 232
        432 41
        85 413
        106 197
        238 73
        96 191
        419 412
        32 390
        134 13
        237 346
        236 345
        437 245
        136 192
        76 404
        35 187
        45 162
        83 104
        93 293
        16 370
        204 175
        9 384
        132 1
        76 275
        242 432
        367 165
        412 318
        130 450
        240 48
        3 88
        429 334
        30 271
        98 215
        141 442
        62 266
        395 5
        12 303
        307 426
        230 245
        159 52
        347 104
        418 434
        209 111
        15 431
        40 43
        339 33
        105 428
        114 311
        16 218
        137 351
        381 397
        246 255
        251 149
        141 258
        74 102
        329 397
        374 131
        162 230
        323 174
        412 280
        66 253
        125 315
        113 420
        273 183
        380 47
        14 416
        321 207
        328 384
        107 211
        142 46
        398 343
        291 15
        250 85
        222 406
        355 256
        454 427
        367 133
        37 286
        142 301
        275 146
        19 309
        159 165
        403 410
        37 167
        181 56
        422 206
        304 333
        267 358
        308 131
        144 197
        36 318
        160 337
        28 310
        68 201
        69 202
        134 321
        141 396
        393 82
        6 316
        255 345
        26 212
        139 77
        317 420
        414 28
        134 104
        428 93
        263 184
        307 246
        400 161
        170 231
        125 172
        360 248
        184 22
        74 361
        108 431
        132 151
        132 72
        251 207
        31 43
        105 85
        57 75
        86 163
        126 17
        267 177
        274 52
        213 123
        26 39
        373 36
        354 355
        112 6
        280 287
        30 166
        138 48
        128 167
        199 404
        173 247
        163 273
        162 42
        158 103
        330 21
        129 315
        155 187
        276 5
        193 440
        113 302
        3 53
        18 424
        181 84
        236 64
        35 227
        388 64
        395 420
        393 421
        286 58
        261 431
        352 10
        112 7
        229 188
        164 261
        404 295
        448 415
        208 169
        298 320
        147 359
        368 205
        121 1
        217 152
        249 68
        411 274
        78 255
        250 69
        60 233
        16 300
        249 66
        188 217
        238 444
        41 148
        428 59
        378 446
        208 382
        149 143
        65 187
        377 365
        168 111
        130 419
        112 193
        287 36
        26 247
        329 182
        347 408
        240 214
        427 159
        375 174
        371 343
        33 296
        208 332
        430 357
        203 398
        227 7
        243 301
        385 365
        379 428
        184 295
        368 195
        19 236
        248 415
        435 24
        279 400
        170 169
        101 234
        203 52
        304 376
        236 319
        237 320
        384 27
        94 160
        444 206
        358 207
        33 308
        254 180
        280 273
        22 2
        164 146
        184 222
        373 352
        78 357
        455 260
        18 286
        287 383
        186 345
        89 358
        153 297
        62 148
        20 168
        439 73
        90 91
        76 23
        381 299
        240 290
        149 455
        292 67
        351 61
        306 109
        133 201
        325 332
        126 161
        189 38
        353 197
        387 306
        387 353
        137 100
        367 392
        90 401
        72 445
        435 64
        124 356
        221 300
        120 323
        423 2
        313 255
        148 270
        426 345
        315 390
        284 358
        1 278
        98 356
        285 359
        283 337
        223 289
        29 196
        411 392
        163 25
        210 341
        277 389
        142 3
        265 179
        14 196
        258 301
        309 281
        372 65
        101 88
        34 150
        125 132
        327 100
        6 120
        324 409
        157 425
        436 63
        251 118
        373 453
        144 72
        406 211
        317 330
        388 186
        250 178
        94 59
        157 374
        34 41
        172 32
        362 216
        19 70
        223 415
        235 286
        348 413
        94 322
        157 290
        269 123
        149 90
        249 312
        274 298
        387 233
        205 223
        242 153
        259 283
        235 49
        44 63
        226 95
        229 156
        120 440
        339 177
        30 333
        137 209
        129 107
        56 299
        154 378
        241 176
        257 168
        106 218
        282 277
        75 295
        417 375
        377 416
        344 225
        228 200
        313 4
        451 317
        8 143
        225 231
        280 145
        369 10
        35 270
        363 99
        146 406
        348 397
        7 65
        332 194
        305 169
        213 147
        357 186
        347 60
        264 276
        344 290
        108 361
        322 202
        126 300
        368 54
        268 196
        150 314
        121 222
        425 350
        333 202
        413 97
        51 407
        139 372
        168 351
        329 362
        77 156
        122 248
        80 294
        320 152
        9 421
        328 164
        291 27
        398 216
        285 48
        83 386
        311 278
        391 402
        281 316
        113 436
        92 224
        226 119
        323 84
        23 50
        391 275
        17 195
        342 9
        56 53
        414 77
        118 8
        408 172
        312 443
        42 217
        112 153
        264 135
        263 185
        395 21
        452 364
        133 292
        202 50
        175 333
        413 299
        342 180
        71 396
        360 296
        431 160
        24 293
        229 411
        407 341
        366 139
        45 248
        121 449
        99 177
        226 453
        436 79
        392 103
        223 409
        452 119
        235 244
        83 272
        117 381
        138 46
        173 123
        131 91
        342 198
        121 82
        319 446
        118 239
        96 87
        265 144
        274 258
        393 278
        196 272
        288 293
        114 226
        190 441
        99 25
        380 419
        230 152
        445 185
        215 224
        140 61
        240 7
        377 272
        0 29
        6 214
        49 331
        379 348
        136 364
        29 239
        281 153
        340 291
        242 147
        415 252
        100 50
        3 323
        243 200
        204 271
        134 179
        113 231
        11 362
        235 67
        214 71
        40 54
        57 306
        268 207
        150 215
        183 4
        64 61
        40 443
        104 325
        385 72
        363 38
        293 174
        32 370
        81 102
        363 310
        448 68
        325 282
        309 39
        322 402
        330 455
        49 194
        51 390
        266 439
        """

        self.waters = """
        0.5 0.89571 0.37061
        0.0 0.62138 0.54386
        0.5 0.51676 0.53777
        0.625 0.26786 0.00495
        0.5 0.60092 0.12529
        0.5 0.12847 0.47296
        0.18426 0.25458 0.1447
        0.18426 0.22435 0.28924
        0.30679 0.97728 0.49448
        0.31821 0.50071 0.43327
        0.68426 0.72435 0.28924
        0.0 0.1575 0.76754
        0.18426 0.94442 0.72588
        0.30926 0.87906 0.55295
        0.69334 0.81517 0.3385
        0.80926 0.37906 0.55295
        0.30926 0.7455 0.72288
        0.5 0.82145 0.76367
        0.81834 0.09977 0.68129
        0.31821 0.40717 0.1071
        0.875 0.58762 0.9859
        0.30679 0.05565 0.40047
        0.5 0.58678 0.57297
        0.19334 0.51528 0.82037
        0.0 0.37153 0.97296
        0.375 0.81469 0.13506
        0.68166 0.40023 0.18129
        0.125 0.41239 0.4859
        0.375 0.99247 0.15246
        0.69334 0.89899 0.41743
        0.68166 0.54436 0.88246
        0.5 0.86201 0.6007
        0.19334 0.65319 0.70362
        0.81834 0.90023 0.18129
        0.80666 0.39899 0.41743
        0.25 0.25165 0.43255
        0.0 0.72406 0.20667
        0.5 0.07077 0.74555
        0.0 0.89935 0.22136
        0.5 0.39935 0.22136
        0.5 0.88188 0.73878
        0.68166 0.34722 0.39955
        0.30679 0.96706 0.09353
        0.5 0.90149 0.66822
        0.69074 0.18634 0.31749
        0.30679 0.8504 0.07509
        0.81574 0.25458 0.1447
        0.68166 0.65278 0.89955
        0.81574 0.22435 0.28924
        0.81574 0.94442 0.72588
        0.125 0.45707 0.8344
        0.0 0.65262 0.82888
        0.69334 0.10102 0.91743
        0.75 0.24836 0.93255
        0.69074 0.91084 0.77646
        0.5 0.30266 0.74918
        0.625 0.23627 0.86189
        0.5 0.6575 0.76754
        0.625 0.00753 0.65246
        0.19074 0.41084 0.77646
        0.0 0.80266 0.74918
        0.0 0.48093 0.99846
        0.75 0.25165 0.43255
        0.625 0.20469 0.40436
        0.0 0.44936 0.05588
        0.30926 0.18634 0.31749
        0.0 0.9035 0.93312
        0.69334 0.01528 0.82037
        0.30679 0.94436 0.90047
        0.80679 0.44436 0.90047
        0.19321 0.3504 0.07509
        0.0 0.16915 0.09328
        0.68179 0.67204 0.51354
        0.18179 0.17204 0.51354
        0.68426 0.32479 0.54093
        0.5 0.60066 0.72136
        0.31834 0.5102 0.75632
        0.5 0.09851 0.16822
        0.0 0.59851 0.16822
        0.81821 0.17204 0.51354
        0.5 0.01908 0.49846
        0.31574 0.32479 0.54093
        0.0 0.51908 0.49846
        0.81574 0.82479 0.54093
        0.25 0.24836 0.93255
        0.68426 0.36963 0.81511
        0.19074 0.7302 0.07614
        0.30666 0.15319 0.70362
        0.68179 0.32796 0.01354
        0.18166 0.81372 0.25331
        0.0 0.93609 0.35495
        0.0 0.95745 0.27158
        0.5 0.45745 0.27158
        0.31821 0.39035 0.89634
        0.0 0.38188 0.73878
        0.19074 0.58916 0.27646
        0.5 0.1319 0.74659
        0.875 0.29531 0.90436
        0.875 0.54294 0.3344
        0.30666 0.84682 0.20362
        0.19321 0.44436 0.90047
        0.80679 0.3504 0.07509
        0.5 0.33085 0.59328
        0.18166 0.15278 0.89955
        0.0 0.83085 0.59328
        0.5 0.36207 0.77065
        0.5 0.7969 0.61426
        0.80666 0.54614 0.709
        0.0 0.2969 0.61426
        0.80666 0.68483 0.8385
        0.81834 0.15278 0.89955
        0.19334 0.60102 0.91743
        0.30926 0.2545 0.22288
        0.81821 0.10965 0.39634
        0.31821 0.60965 0.39634
        0.69321 0.03295 0.59353
        0.68426 0.63037 0.31511
        0.68426 0.24542 0.6447
        0.375 0.91239 0.4859
        0.5 0.63794 0.27065
        0.30926 0.2302 0.07614
        0.0 0.55064 0.55588
        0.5 0.94936 0.05588
        0.0 0.3681 0.24659
        0.68179 0.60965 0.39634
        0.875 0.68532 0.63506
        0.30926 0.81367 0.81749
        0.69334 0.98472 0.32037
        0.18166 0.09977 0.68129
        0.68166 0.59977 0.68129
        0.80926 0.62095 0.05295
        0.81834 0.98981 0.25632
        0.80679 0.6496 0.57509
        0.30666 0.01528 0.82037
        0.18426 0.82479 0.54093
        0.81821 0.09283 0.6071
        0.80926 0.58916 0.27646
        0.31821 0.4993 0.93327
        0.69074 0.2545 0.22288
        0.5 0.11812 0.23878
        0.80679 0.52272 0.99448
        0.81574 0.07375 0.09919
        0.69074 0.2302 0.07614
        0.18179 0.00071 0.43327
        0.625 0.73214 0.50495
        0.625 0.81469 0.13506
        0.19074 0.44263 0.6621
        0.68166 0.31372 0.25331
        0.625 0.30636 0.45599
        0.18166 0.95564 0.38246
        0.68166 0.45564 0.38246
        0.68179 0.59283 0.6071
        0.0 0.98324 0.03777
        0.31834 0.31372 0.25331
        0.68179 0.4993 0.93327
        0.30679 0.15697 0.44554
        0.5 0.13799 0.1007
        0.0 0.08901 0.21071
        0.0 0.13295 0.85733
        0.81834 0.04436 0.88246
        0.0 0.40149 0.66822
        0.375 0.79531 0.90436
        0.18179 0.90717 0.1071
        0.31574 0.75458 0.1447
        0.31574 0.42626 0.59919
        0.0 0.06392 0.85495
        0.5 0.56392 0.85495
        0.30666 0.04614 0.709
        0.125 0.58762 0.9859
        0.0 0.01676 0.53777
        0.0 0.02622 0.45656
        0.5 0.52622 0.45656
        0.125 0.68532 0.63506
        0.0 0.42923 0.24555
        0.31821 0.32796 0.01354
        0.0 0.5521 0.81313
        0.31574 0.67521 0.04093
        0.5 0.8681 0.24659
        0.5 0.4035 0.93312
        0.30926 0.76981 0.57614
        0.19334 0.48472 0.32037
        0.375 0.23627 0.86189
        0.31574 0.27565 0.78924
        0.5 0.66915 0.09328
        0.31821 0.59283 0.6071
        0.31821 0.67204 0.51354
        0.875 0.49247 0.15246
        0.375 0.20469 0.40436
        0.18426 0.07375 0.09919
        0.0 0.8425 0.26754
        0.19334 0.39899 0.41743
        0.69334 0.15319 0.70362
        0.0 0.61812 0.23878
        0.5 0.22406 0.20667
        0.69074 0.94263 0.6621
        0.69074 0.81367 0.81749
        0.81834 0.84722 0.39955
        0.69074 0.76981 0.57614
        0.5 0.43609 0.35495
        0.0 0.57077 0.74555
        0.125 0.19365 0.95599
        0.375 0.95707 0.8344
        0.875 0.45707 0.8344
        0.5 0.10429 0.87061
        0.0 0.60429 0.87061
        0.81821 0.89035 0.89634
        0.30679 0.03295 0.59353
        0.18166 0.84722 0.39955
        0.81821 0.9851 0.55543
        0.31834 0.54436 0.88246
        0.5 0.63295 0.85733
        0.875 0.50753 0.65246
        0.68179 0.40717 0.1071
        0.80666 0.34682 0.20362
        0.0 0.21665 0.15724
        0.80666 0.48472 0.32037
        0.30666 0.18483 0.8385
        0.18179 0.01491 0.05543
        0.5 0.77594 0.70667
        0.5 0.92923 0.24555
        0.625 0.18532 0.63506
        0.125 0.73627 0.86189
        0.19321 0.53295 0.59353
        0.69321 0.84304 0.94554
        0.68166 0.48981 0.25632
        0.18179 0.10965 0.39634
        0.31574 0.63037 0.31511
        0.125 0.26374 0.36189
        0.18426 0.17521 0.04093
        0.30926 0.12095 0.05295
        0.0 0.91323 0.07297
        0.0 0.09651 0.43312
        0.5 0.59651 0.43312
        0.81574 0.77565 0.78924
        0.875 0.31469 0.13506
        0.81834 0.0102 0.75632
        0.19321 0.46706 0.09353
        0.69321 0.96706 0.09353
        0.30679 0.1496 0.57509
        0.625 0.91239 0.4859
        0.0 0.19734 0.24918
        0.375 0.69365 0.95599
        0.5 0.3425 0.26754
        0.875 0.19365 0.95599
        0.0 0.04255 0.77158
        0.0 0.86703 0.0079
        0.31574 0.55559 0.22588
        0.80666 0.45387 0.209
        0.5 0.87862 0.04386
        0.18179 0.89035 0.89634
        0.68179 0.39035 0.89634
        0.30666 0.89899 0.41743
        0.30679 0.84304 0.94554
        0.0 0.97379 0.95656
        0.125 0.54294 0.3344
        0.19074 0.55737 0.1621
        0.19074 0.26981 0.57614
        0.19074 0.62095 0.05295
        0.69074 0.12095 0.05295
        0.5 0.39908 0.62529
        0.5 0.9479 0.31313
        0.19074 0.37906 0.55295
        0.0 0.4479 0.31313
        0.19321 0.6496 0.57509
        0.69321 0.1496 0.57509
        0.375 0.73214 0.50495
        0.875 0.23214 0.50495
        0.5 0.84738 0.32888
        0.0 0.86705 0.35733
        0.0 0.34738 0.32888
        0.375 0.30636 0.45599
        0.80666 0.60102 0.91743
        0.875 0.80636 0.45599
        0.5 0.71665 0.15724
        0.625 0.08762 0.9859
        0.31574 0.44442 0.72588
        0.5 0.12138 0.54386
        0.30926 0.94263 0.6621
        0.0 0.62847 0.47296
        0.125 0.76786 0.00495
        0.68426 0.75458 0.1447
        0.19334 0.34682 0.20362
        0.18426 0.92626 0.59919
        0.68426 0.42626 0.59919
        0.375 0.76374 0.36189
        0.875 0.26374 0.36189
        0.69334 0.04614 0.709
        0.80926 0.7545 0.22288
        0.125 0.29531 0.90436
        0.625 0.79531 0.90436
        0.0 0.13794 0.27065
        0.875 0.41239 0.4859
        0.5 0.0521 0.81313
        0.19321 0.34304 0.94554
        0.5 0.05064 0.55588
        0.31834 0.59977 0.68129
        0.81821 0.90717 0.1071
        0.19334 0.31517 0.3385
        0.69321 0.02272 0.99448
        0.68426 0.27565 0.78924
        0.18426 0.77565 0.78924
        0.81574 0.17521 0.04093
        0.81574 0.13037 0.31511
        0.18166 0.0102 0.75632
        0.68166 0.5102 0.75632
        0.0 0.08678 0.57297
        0.68166 0.68628 0.75331
        0.31834 0.48981 0.25632
        0.69334 0.95387 0.209
        0.31834 0.40023 0.18129
        0.30666 0.95387 0.209
        0.19321 0.65697 0.44554
        0.18426 0.86963 0.81511
        0.31574 0.57375 0.09919
        0.68179 0.50071 0.43327
        0.80666 0.65319 0.70362
        0.125 0.31469 0.13506
        0.625 0.04294 0.3344
        0.0 0.70311 0.11426
        0.31821 0.51491 0.05543
        0.81821 0.01491 0.05543
        0.125 0.80636 0.45599
        0.80926 0.41084 0.77646
        0.375 0.26786 0.00495
        0.875 0.76786 0.00495
        0.0 0.89908 0.62529
        0.18179 0.9993 0.93327
        0.0 0.43611 0.94403
        0.31821 0.4851 0.55543
        0.19074 0.2455 0.72288
        0.375 0.04294 0.3344
        0.0 0.91099 0.71071
        0.81574 0.92626 0.59919
        0.80666 0.51528 0.82037
        0.31574 0.72435 0.28924
        0.0 0.39571 0.37061
        0.625 0.95707 0.8344
        0.80926 0.44263 0.6621
        0.68179 0.51491 0.05543
        0.69334 0.84682 0.20362
        0.80679 0.47728 0.49448
        0.31834 0.65278 0.89955
        0.31834 0.45564 0.38246
        0.69334 0.18483 0.8385
        0.18426 0.13037 0.31511
        0.125 0.49247 0.15246
        0.625 0.99247 0.15246
        0.0 0.78336 0.65724
        0.0 0.32145 0.76367
        0.5 0.28336 0.65724
        0.18166 0.98981 0.25632
        0.19321 0.52272 0.99448
        0.80926 0.68634 0.31749
        0.81574 0.74542 0.6447
        0.375 0.18532 0.63506
        0.31574 0.24542 0.6447
        0.80679 0.55565 0.40047
        0.80926 0.55737 0.1621
        0.30666 0.81517 0.3385
        0.80666 0.31517 0.3385
        0.69321 0.8504 0.07509
        0.80926 0.26981 0.57614
        0.18166 0.18628 0.75331
        0.18166 0.90023 0.18129
        0.68426 0.55559 0.22588
        0.875 0.70469 0.40436
        0.30926 0.08916 0.27646
        0.18166 0.04436 0.88246
        0.81574 0.86963 0.81511
        0.5 0.69734 0.24918
        0.31834 0.68628 0.75331
        0.81834 0.18628 0.75331
        0.5 0.17855 0.26367
        0.0 0.67855 0.26367
        0.81574 0.05559 0.22588
        0.5 0.36703 0.0079
        0.5 0.54255 0.77158
        0.75 0.75165 0.43255
        0.5 0.47379 0.95656
        0.19074 0.31367 0.81749
        0.625 0.69365 0.95599
        0.80926 0.2455 0.72288
        0.69321 0.97728 0.49448
        0.81834 0.81372 0.25331
        0.19321 0.47728 0.49448
        0.80679 0.65697 0.44554
        0.69074 0.87906 0.55295
        0.69074 0.7455 0.72288
        0.80679 0.46706 0.09353
        0.375 0.00753 0.65246
        0.0 0.6319 0.74659
        0.5 0.41099 0.71071
        0.30666 0.10102 0.91743
        0.0 0.56389 0.44403
        0.0 0.86207 0.77065
        0.5 0.06389 0.44403
        0.0 0.10092 0.12529
        0.0 0.27594 0.70667
        0.5 0.15262 0.82888
        0.81821 0.00071 0.43327
        0.25 0.74836 0.93255
        0.81834 0.95564 0.38246
        0.68426 0.44442 0.72588
        0.25 0.75165 0.43255
        0.19334 0.54614 0.709
        0.5 0.98093 0.99846
        0.125 0.50753 0.65246
        0.19334 0.68483 0.8385
        0.18426 0.74542 0.6447
        0.81821 0.82796 0.01354
        0.125 0.70469 0.40436
        0.375 0.08762 0.9859
        0.80926 0.7302 0.07614
        0.80926 0.31367 0.81749
        0.30926 0.05737 0.1621
        0.5 0.87153 0.97296
        0.625 0.76374 0.36189
        0.5 0.41323 0.07297
        0.875 0.73627 0.86189
        0.68426 0.67521 0.04093
        0.69321 0.05565 0.40047
        0.19321 0.55565 0.40047
        0.18179 0.9851 0.55543
        0.68179 0.4851 0.55543
        0.0 0.10066 0.72136
        0.18426 0.05559 0.22588
        0.19334 0.45387 0.209
        0.81821 0.9993 0.93327
        0.31574 0.36963 0.81511
        0.19074 0.7545 0.22288
        0.68426 0.57375 0.09919
        0.0 0.36201 0.6007
        0.5 0.36705 0.35733
        0.30679 0.02272 0.99448
        0.75 0.74836 0.93255
        0.0 0.37862 0.04386
        0.69321 0.15697 0.44554
        0.18179 0.82796 0.01354
        0.80679 0.34304 0.94554
        0.125 0.23214 0.50495
        0.5 0.20311 0.11426
        0.31834 0.34722 0.39955
        0.69074 0.05737 0.1621
        0.30926 0.91084 0.77646
        0.18179 0.09283 0.6071
        0.5 0.63297 0.5079
        0.5 0.48324 0.03777
        0.0 0.13297 0.5079
        0.5 0.93611 0.94403
        0.80679 0.53295 0.59353
        0.0 0.63799 0.1007
        0.69074 0.08916 0.27646
        0.5 0.58901 0.21071
        0.19074 0.68634 0.31749
        0.69321 0.94436 0.90047
        0.30666 0.98472 0.32037
        """

        self.coord = "relative"

        self.cages = """
        14 0.23703 0.15911 0.17915
        12 -0.5 -0.05041 0.57473
        15 1.0 1.04328 0.32986
        12 0.77337 0.89599 0.30088
        16 0.5 0.40908 0.48818
        12 0.22663 0.10401 0.80088
        12 -0.27285 -0.42375 0.5041
        12 -0.5 -0.02541 0.73983
        12 0.22715 0.92375 1.0041
        12 0.5 0.02579 0.91092
        14 0.5 0.30976 0.89442
        16 -0.5 -0.31086 0.62899
        12 1.0 0.52541 0.23983
        16 1.0 0.81086 0.12899
        15 -0.5 -0.54328 0.82986
        12 0.22663 0.89599 1.30088
        14 -0.23703 -0.15911 0.67915
        12 0.0 0.44959 0.57473
        12 0.77337 0.10401 0.80088
        12 0.27285 0.42375 0.0041
        14 0.5 0.69024 0.39442
        15 0.0 -0.31666 0.96975
        14 0.23703 -0.15911 0.67915
        14 0.0 0.19024 0.39442
        12 0.5 0.97421 0.41092
        15 -0.5 -0.49454 0.65469
        12 1.0 0.60048 0.6406
        12 0.5 0.05041 0.07473
        12 0.0 0.37116 0.85699
        12 -0.27285 0.42375 1.0041
        15 0.0 -0.26747 0.52665
        12 -0.27337 0.39599 1.30088
        15 0.5 0.76747 0.02665
        14 1.0 0.76715 0.33939
        12 0.77285 0.92375 0.0041
        14 -0.5 -0.26715 0.83939
        12 1.0 0.62884 0.35699
        15 0.5 0.49454 0.15469
        12 -0.5 -0.10048 1.1406
        14 0.73703 0.65911 0.17915
        14 0.5 0.26715 0.33939
        12 0.77285 0.07625 0.5041
        15 0.5 0.54328 0.32986
        16 0.0 0.09092 0.98818
        12 0.27337 -0.39599 0.80088
        14 0.73703 0.34089 0.67915
        15 0.0 0.31666 0.46975
        12 0.5 0.02541 0.23983
        12 0.5 0.10048 0.6406
        12 0.0 0.47459 0.73983
        12 -0.5 -0.12884 0.85699
        12 0.0 -0.47421 0.91092
        12 0.22715 0.07625 0.5041
        12 0.0 -0.28089 0.74399
        16 0.0 0.18914 0.62899
        14 0.26297 0.65911 1.17915
        14 0.0 0.23285 0.83939
        12 0.0 0.39952 1.1406
        16 -0.5 -0.40908 0.98818
        12 0.27285 -0.42375 0.5041
        15 0.0 -0.04328 0.82986
        15 0.5 0.18334 0.96975
        16 1.0 0.90908 0.48818
        15 1.0 0.99454 0.15469
        12 0.27337 0.39599 0.30088
        15 0.5 0.81666 0.46975
        12 0.0 0.47421 0.41092
        15 0.0 0.00546 0.65469
        15 0.5 0.23253 0.52665
        14 0.26297 0.34089 0.67915
        12 1.0 0.55041 0.07473
        14 0.0 -0.19024 0.89442
        12 -0.27337 -0.39599 0.80088
        12 0.5 0.21911 0.74399
        14 -0.23703 0.15911 1.17915
        12 0.0 0.28089 0.24399
        16 0.5 0.31086 0.12899
        12 0.5 0.12884 0.35699
        12 0.5 0.78089 0.24399
        15 0.0 0.26747 0.02665
        """

        self.bondlen = 3

        self.cell = """
        14.589965896299718 50.38172056030764 39.28872908930472
        """

        self.density = 0.4719558184140425

        self.cell = cellvectors(a=14.589965896299718,
                                b=50.38172056030764,
                                c=39.28872908930472)
